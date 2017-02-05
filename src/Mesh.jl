# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: string, print, show

# Export types
export Mesh

# Export methods
export delete, string, print, show, load, save,
  numnodes, numelems, numdims, boundingbox, fullsize

# Type definitions
type Mesh

  ptr::Cxx.CppPtr
  rowptr::Array{Cint,1}
  colidx::Array{Cint,1}

  function Mesh(ptr::Cxx.CppPtr)

    # Perform initialisation
    _info("Initiliasing mesh.")
    @cxx ptr->Setup()

    # Compute sparsity pattern
    _info("Calculating sparsity pattern.")
    rowptr, colidx = _sparsity!(ptr)

    mesh = new(ptr, rowptr, colidx)
    finalizer(mesh, _mesh_delete)

    return mesh

  end

end

"""
    Mesh(node, elem, eltp)

Construct a new mesh given a list of vertices, elements, and element types.
"""
function Mesh(node::Matrix{Float64}, elem::Matrix{Float64}, eltp::Vector{Float64})
  meshptr = _mesh_new()
  _make!(meshptr,node,elem,eltp)
  return Mesh(meshptr)
end

"""
    Mesh(filename)

Construct a new mesh from a specified file.
"""
function Mesh(fn::String)
  meshptr = _mesh_new()
  _load!(meshptr,fn)
  return Mesh(meshptr)
end

# Create new mesh pointer
_mesh_new() = @cxxnew TOAST_Mesh()

# Delete and deallocate mesh pointer
# TODO: Does Cxx.jl do this for me?
_mesh_delete(mesh) = finalize(mesh.ptr) # icxx"""delete $(mesh.ptr);"""

# Compute sparsity structure in zero-based CSR format and 1-based CSC.
function _sparsity!(ptr::Cxx.CppPtr)

  rowptr = Ptr{Int32}[0]
  colidx = Ptr{Int32}[0]
  nnzero = Int32[0]

  icxx"""
    $(ptr)->SparseRowStructure($(pointer(rowptr))[0], $(pointer(colidx))[0], $(pointer(nnzero))[0]);
  """

  nnd = @cxx (@cxx ptr->nlist)->Len()

  _info("Number of vertices: ", nnd, ", number of nonzeros: ", nnzero[1])

  # Zero based rowptr and colidx, for resuse when building a system matrix
  rowptr = unsafe_wrap(Array,rowptr[1],nnd+1,true)
  colidx = unsafe_wrap(Array,colidx[1],nnzero[1],true)

  return rowptr, colidx

end

# Build a mesh from a matrix of vertices, elements, and element types
function _make!(meshptr::Cxx.CppPtr, node, elem, eltp)

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

      Mesh *mesh = $(meshptr);

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

  return nothing

end

# Method definitions
# data
# surf
# maxnodes

# String conversion
string(mesh::Mesh) = ("")

# Print method
print(io::IO, mesh::Mesh) = print(io, string(mesh))

# Show method
show(io::IO, mesh::Mesh) = print(io, mesh)

# Load a mesh from the specified filename
function _load!(ptr::Cxx.CppPtr, fn::String)

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
    ifs >> *($(ptr));
    ifs.close();
  """

  _info("Loaded ", stat(fn).size, " bytes from ", fn, ".")

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

  _info("Mesh written to ", fn, ", ", stat(fn).size, "bytes.")

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
  bb = Array{Cdouble}(dim,2)

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
