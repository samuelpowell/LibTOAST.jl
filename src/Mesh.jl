# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports and exports
import Base: string, print, show

# Export types
export Mesh, SystemMatrix

# Export enumerations
export FF, DD, BNDFF, PFF, PDD, BNDPFF, P, PF, BNDPF

# Export methods
export delete, string, print, show, sparsity!, load, load!, save,
  numnodes, numelemes, numdims, boundingbox, fillsize, make!, assemble!, sysmat

# Type definitions
type Mesh

  # Pointer to C++ mesh object
  ptr::Cxx.CppPtr

  # Initialisation
  initialised::Bool

  # Sparsity pattern
  sparsity
  sparsity_direct
  sparsity_rowptr::Array{Int32,1}
  sparsity_colidx::Array{Int32,1}

  # Constructor
  function Mesh()
    mesh = new((@cxxnew TOAST_Mesh()), false);
    finalizer(mesh, delete)
    return mesh
  end

end


type SystemMatrix

  ptr::Cxx.CppPtr

  function SystemMatrix(mesh)
    sysmat = new(@cxxnew RCompRowMatrix(numnodes(mesh),numnodes(mesh),
      pointer(mesh.sparsity_rowptr), pointer(mesh.sparsity_colidx)))
    finalizer(sysmat, delete)
    return sysmat
  end

end


"""
    Mesh(node, elem, eltp)

Construct a new mesh given a list of vertices, elements, and element types.
"""
function Mesh(node::Matrix{Float64}, elem::Matrix{Float64}, eltp::Vector{Float64})
  mesh = Mesh()
  make!(mesh,node,elem,eltp)
  return mesh
end

"""
    Mesh(filename)

Construct a new mesh from a specified file.
"""
function Mesh(fn::String)
  mesh = Mesh()
  load!(mesh,fn)
  return mesh
end

# Method definitions
# data
# surf
# maxnodes

# Delete the underlying pointer to a Toast++ mesh
delete(mesh::Mesh) = icxx"""delete $(mesh.ptr);"""

# Delete underlying pointer to a Toast++ system matrix
delete(sysmat::SystemMatrix) = icxx"""delete $(sysmat.ptr);"""

# String conversion
string(mesh::Mesh) = ("")

# Print method
print(io::IO, mesh::Mesh) = print(io, string(mesh))

# Show method
show(io::IO, mesh::Mesh) = print(io, mesh)

"""
    sparsity!(mesh)

Compute and sparsity structure in both Toast++ (CSR) and Julia (CSC) formats
and store in mesh object.
"""
function sparsity!(mesh::Mesh)

  rowptr = Ptr{Int32}[0]
  colidx = Ptr{Int32}[0]
  nnzero = Int32[0]

  info("Computing sparsity pattern.")

  icxx"""
    $(mesh.ptr)->SparseRowStructure($(pointer(rowptr))[0], $(pointer(colidx))[0], $(pointer(nnzero))[0]);
  """

  nnd = numnodes(mesh)

  info("Number of nodes: ", nnd, ", Number of nonzeros: ", nnzero[1])

  # Zero based rowptr and colidx, stored for resuse when building, e.g, the
  # system matrix
  mesh.sparsity_rowptr = pointer_to_array(rowptr[1],nnd+1,true)
  mesh.sparsity_colidx = pointer_to_array(colidx[1],nnzero[1],true)

  # Store 1-based Julia CSC sparsity pattern
  mesh.sparsity = _CSC(nnd,nnd,mesh.sparsity_rowptr,mesh.sparsity_colidx,ones(nnzero[1]))

  return nothing

end

"""
    load(filename)

Construct and initialise a Toast++ mesh from a specified file.
"""
load(fn::String) = load!(Mesh(), fn)

"""
    load!(mesh, filename)

Construct and initialise a Toast++ mesh from a specified file.
"""
function load!(mesh::Mesh, fn::String)

  # TODO: We would like to do something like this:
  # ifs = @cxx std::ifstream(pointer(filename))
  # icxx"""$(ifs) >> *$(obj.meshptr);"""

  # Check if file exists
  if !isfile(fn)
    error("File not found ", fn)
  end

  # Load the mesh
  icxx"""
    ifstream ifs($(pointer(fn)));
    ifs >> *($(mesh.ptr));
    ifs.close();
  """

  info("Loaded ", stat(fn).size, " bytes from ", fn, ".")

  # Perform initialisation
  info("Initiliasing mesh.")
  @cxx mesh.ptr->Setup()

  # Calcualte the sparsity pattern
  info("Calculating sparsity pattern.")
  sparsity!(mesh)

  mesh.initialised = true

  return

end

"""
    save(mesh, filename)

Save an initialised mesh to the specified filename.
"""
function save(mesh::Mesh, fn::String)

  # Output the mesh
  icxx"""
    ofstream ofs($(pointer(fn)));
    ofs << *($(mesh.ptr));
    ofs.close();
  """

  info("Mesh written to ", fn, ", ", stat(fn).size, "bytes.")

  return nothing

end

"""
    numnodes(mesh)

Return the number of nodes in the mesh.
"""
numnodes(mesh::Mesh) = @cxx (@cxx mesh.ptr->nlist)->Len()

"""
    numelems(mesh)

Return the number of elements in the mesh.
"""
numelems(mesh::Mesh) = @cxx (@cxx mesh.ptr->elist)->Len()

"""
    numdims(mesh)

Return the number of spatial dimensions of the mesh.
"""
numdims(mesh::Mesh) = @cxx mesh.ptr->Dimension()

"""
    boundingbox(mesh)

Return the bounding box of the mesh in a ``d \times 2`` matrix.
"""
function boundingbox(mesh::Mesh)

  dim = numdims(mesh)
  bb = Array{Float64}(dim,2)

  icxx"""
    double *jptr = $(pointer(bb));

    Point pmin($dim), pmax($dim);
    $(mesh.ptr)->BoundingBox(pmin, pmax);

    for(int i=0; i<$dim; i++) {
      *jptr++ = pmin[i];
      *jptr++ = pmax[i];
    }
  """

  return bb

end

"""
    size(mesh)

Compute the total area (2D) for volume (3D) of the mesh.
"""
fullsize(mesh::Mesh) = @cxx mesh.ptr->FullSize()

"""
    make!(mesh, node, elem, eltp)

Construct a mesh given a matrix of vertices, elements, and a vector of element
types.
"""
function make!(mesh::Mesh, node, elem, eltp)

  jnvtx = size(node,1)
  jnel = size(elem,1)
  jdim = size(node,2)
  jnnd0 = size(elem,2)

  # The following code fragment is modified from mtMesh.cc
  icxx"""
      int i, j, k;
      int nvtx = $jnvtx;
      int nel  = $jnel;
      int dim  = $jdim;
      int nnd0 = $jnnd0;
      double *vtx = $(pointer(node));
      double *idx = $(pointer(elem));
      double *etp = $(pointer(eltp));

      Mesh *mesh = $(mesh.ptr);

      // Create node list
      mesh->nlist.New (nvtx);
      for (i = 0; i < nvtx; i++) {
        mesh->nlist[i].New(dim);
        mesh->nlist[i].SetBndTp (BND_NONE); // Unknown (yet)
      }
      for (j = k = 0; j < dim; j++) {
        for (i = 0; i < nvtx; i++) {
          mesh->nlist[i][j] = vtx[k++];
        }
      }

      // Create element list
      Element *el, **list = new Element*[nel];
      for (i = 0; i < nel; i++) {
        int eltp = (int)(etp[i]+0.5);
        switch (eltp) {
          case ELID_TRI3OLD:
            list[i] = new Triangle3old;
            break;
          case ELID_TET4:
            list[i] = new Tetrahedron4;
            break;
          case ELID_WDG6:
            list[i] = new Wedge6;
            break;
          case ELID_VOX8:
            list[i] = new Voxel8;
            break;
          case ELID_TRI6:
            list[i] = new Triangle6;
            break;
          case ELID_TET10:
            list[i] = new Tetrahedron10;
            break;
          case ELID_TRI6_IP:
            list[i] = new Triangle6_ip;
            break;
          case ELID_TRI10:
            list[i] = new Triangle10;
            break;
          case ELID_TRI10_IP:
            list[i] = new Triangle10_ip;
            break;
          case ELID_TET10_IP:
            list[i] = new Tetrahedron10_ip;
            break;
          case ELID_PIX4:
            list[i] = new Pixel4;
            break;
          case ELID_TRI3:
            list[i] = new Triangle3;
            break;
          case ELID_TRI3D3:
            list[i] = new Triangle3D3;
            break;
          case ELID_TRI3D6:
            list[i] = new Triangle3D6;
            break;
          default:
            printf("Element type not supported!\n");
            list[i] = 0;
            break;
        }
      }
      mesh->elist.SetList (nel, list);
      delete []list;

      for (j = k = 0; j < nnd0; j++) {
        for (i = 0; i < nel; i++) {
          if ((el = mesh->elist[i])) {
            if (j < el->nNode())
              el->Node[j] = (int)(idx[k]-0.5);
            }
          k++;
        }
      }

      // Create dummy parameter list
      mesh->plist.New (nvtx);
      mesh->plist.SetMua (0.01);
      mesh->plist.SetMus (1);
      mesh->plist.SetN (1);

      // Setup mesh
      mesh->Setup();
    """

  sparsity!(obj)

  obj.initialised = true

  print(obj)

  return nothing

end

# The following enumerations are defined in mesh.h
@enum BilinearIntegrals FF=0 DD=1 BNDFF=12
@enum BilinearParamIntegrals PFF=2 PDD=3 BNDPFF=4
@enum LinearParamIntegrals P=0 PF=1 BNDPF=2

"""
    assemble!(sysmat, mesh, integral)

Assemble the bilinear form over the mesh as specified by integral, and add the
result to the provided system matrix.
"""
function assemble!(sysmat::SystemMatrix,
                   mesh::Mesh,
                   integral::BilinearIntegrals)

  icxx"""AddToSysMatrix (*$(mesh.ptr), *$(sysmat.ptr), NULL, $(Cint(integral)));"""

  return nothing

end

"""
    assemble!(sysmat, mesh, integral, param)

Assemble the bilinear form over the mesh as specified by integral and associated
parameter, add the result to the provided system matrix.
"""
function assemble!(sysmat::SystemMatrix,
                   mesh::Mesh,
                   integral::BilinearParamIntegrals,
                   param::Vector{Float64})

  nprm = length(param)
  pprm = pointer(param)
  mode = Cint(integral)

  if nprm != numnodes(mesh)
    error("Parameter vector length ($nprm) not equal to nodal basis ($(numnodes(mesh))).")
  end

  icxx"""
    RVector prm($nprm, $pprm, SHALLOW_COPY);
    AddToSysMatrix (*$(mesh.ptr), *$(sysmat.ptr), &prm, $mode);
  """

  return nothing

end

function sysmat(F::SystemMatrix, mesh::Mesh)
  nnd = numnodes(mesh)
  _CSC(nnd,nnd,mesh.sparsity_rowptr,mesh.sparsity_colidx, @cxx F.ptr->ValPtr())
end
