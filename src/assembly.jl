# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Export integral enumerations (see mesh.h)
export BilinearIntegrals, FF, DD, BNDFF
export BilinearParamIntegrals, PFF, PDD, BNDPFF
export LinearParamIntegrals, P, PF, BNDPF

# Export special purpose assembly types
export PointSource

# Methods
export assemble, assemble!

# The following enumerations are defined in mesh.h
"""
  BilinearIntegrals

An enumeration of integrals over products of two basis functions, supported for
sytsem matrix assembly.

# Integrals
* `FF`: ``∫ uᵢ(r) uⱼ(r) dΩ``
* `DD`: ``∫ ∇uᵢ(r) ⋅ ∇uⱼ(r) dΩ``
* `BNDFF`: ``∫ uᵢ(r) uⱼ(r) d(δΩ)``

Where ``uᵢ(r)`` is the basis function for the ``i``th node, and ``d(δΩ)``
indicates integration only over boundary nodes.
"""
@enum BilinearIntegrals FF=0 DD=1 BNDFF=12

"""
  BilinearParamIntegrals

An enumeration of integrals over products of two basis functions and a function
defined in the nodal basis, supported for sytsem matrix assembly.

# Integrals
* `PFF`: ``∑ₖ fₖ ∫ uᵢ(r) uⱼ(r) uₖ(r) dΩ``
* `PDD`: ``∑ₖ fₖ ∫ ∇uᵢ(r) ⋅ ∇uⱼ(r) uₖ(r) dΩ``
* `BNPDFF`: ``∑ₖ fₖ ∫ uᵢ(r) uⱼ(r) uₖ(r) d(δΩ)``

Where ``uᵢ(r)`` is the basis function for the ``i``th node, and ``d(δΩ)``
indicates integration only over boundary nodes.
"""
@enum BilinearParamIntegrals PFF=2 PDD=3 BNDPFF=4

"""
  BilinearParamIntegrals

An enumeration of integrals over a function defined in the nodal basis, optionally
in product with a basis function, supported for RHS assembly.

# Integrals
* `P`: ``∑ₖ fₖ ∫ uₖ(r) dΩ``
* `PF`: ``∑ₖ fₖ ∫ uᵢ(r) uₖ(r) dΩ``
* `BNPDF`: ``∑ₖ fₖ ∫ uᵢ(r) uₖ(r) d(δΩ)``

Where ``uᵢ(r)`` is the basis function for the ``i``th node, and ``d(δΩ)``
indicates integration only over boundary nodes.
"""
@enum LinearParamIntegrals P=0 PF=1 BNDPF=2

"""
  PointSource

A singleton type used to dispatch the assembly of a point source right hand side
"""
type PointSource end

"""
    assemble(mesh, integral)

Construct a new system matrix and assemble the bilinear form over the mesh,
as specified by integral.

# Arguments
* `mesh::Mesh`: a TOAST mesh
* `integral::BilinearIntegrals`: the desired integral
"""
function assemble(mesh::Mesh, int::BilinearIntegrals)
  sysmat = SystemMatrix(mesh)
  assemble!(sysmat, int)
  return sysmat
end

"""
    assemble!(sysmat, integral)

Assemble the bilinear form over the mesh as specified by integral, and add the
result to the provided system matrix.

# Arguments
* `sysmat::SystemMatrix`: a TOAST system matrix
* `integral::BilinearIntegrals`: the desired integral
"""
function assemble!(sysmat::SystemMatrix, int::BilinearIntegrals)
  meshptr = sysmat.mesh.ptr
  icxx"""AddToSysMatrix (*$(meshptr), *$(sysmat.ptr), NULL, $(Cint(int)));"""
  return nothing
end

"""
    assemble(mesh, integral, parameter)

Construct a new system matrix and assemble the bilinear form over the mesh,
as specified by integral and associated parameter.

# Arguments
* `mesh::Mesh`: a TOAST mesh
* `integral::BilinearIntegrals`: the desired integral
* `parameter::NodalCoeff`: a function in the nodal basis
"""
function assemble(mesh::Mesh, int::BilinearParamIntegrals, param::NodalCoeff)
  sysmat = SystemMatrix(mesh)
  assemble!(sysmat, int, param)
  return sysmat
end

"""
    assemble!(sysmat, integral, parameter)

Assemble the bilinear form over the mesh as specified by integral and associated
parameter, add the result to the provided system matrix.

# Arguments
* `sysmat::SystemMatrix`: a TOAST system matrix
* `integral::BilinearIntegrals`: the desired integral
* `parameter::NodalCoeff`: a function in the nodal basis
"""
function assemble!(sysmat::SystemMatrix,
                   int::BilinearParamIntegrals,
                   param::NodalCoeff)

  assert(param.mesh == sysmat.mesh)

  nprm = length(param)
  pprm = pointer(param.data)
  mode = Cint(int)
  nnd = nodecount(sysmat.mesh)

  if nprm != nnd
    error("Parameter vector length ($nprm) not equal to nodal basis ($nnd).")
  end

  icxx"""
    RVector prm($nprm, $pprm, SHALLOW_COPY);
    AddToSysMatrix (*$(sysmat.mesh.ptr), *$(sysmat.ptr), &prm, $mode);
  """

  return nothing

end

"""
    assemble(mesh, integral, parameter)

Construct a new RHS and assemble the linear form over the `mesh`, as specified
by the `integral` and associated `parameter`.

# Arguments
* `mesh::Mesh`: a TOAST mesh
* `integral::LinearParamIntegrals`: the desired integral
* `parameter::NodalCoeff`: a function in the nodal basis
"""
function assemble(mesh::Mesh, int::LinearParamIntegrals, param::NodalCoeff)
  rhs = zero(NodalCoeff, mesh)
  assemble!(rhs, int, param)
  return rhs
end

"""
    assemble!(rhs, integral, parameter)

Assemble the bilinear form over the mesh as specified by the `integral` and
associated `parameter`, add the result to the provided NodalCoeff vector `rhs`.

# Arguments
* `rhs::NodalCoeff`: a nodal coefficient vector
* `integral::LinearParamIntegrals`: the desired integral
* `parameter::NodalCoeff`: a function in the nodal basis
"""
function assemble!(rhs::NodalCoeff,
                   int::LinearParamIntegrals,
                   param::NodalCoeff)

  assert(param.mesh == rhs.mesh)

  nprm = length(param)
  pprm = pointer(param.data)
  prhs = pointer(rhs.data)
  mode = Cint(int)
  nnd = nodecount(rhs.mesh)

  if nprm != nnd
    error("Parameter vector length ($nprm) not equal to nodal basis ($nnd).")
  end

  #
  # icxx"""
  #   RVector prm($nprm, $pprm, SHALLOW_COPY);
  #   RVector rhs($nprm, $prhs, SHALLOW_COPY);
  #   AddToRHS_elasticity (*$(rhs.mesh.ptr), rhs, &prm, $mode);
  # """

  error("Linear integral assembly not implemented")

  return nothing

end

"""
    assemble(mesh, PointSource, centre)

Special purpose assembly function to generate an right hand side source vector
representing a unitary isotropic point source at the specified centre location.
"""
function assemble(mesh, ::Type{PointSource}, centre)

  rhs = zero(NodalCoeff, mesh)

  prhs = pointer(rhs.data)
  nvtx = nodecount(mesh)
  pcnt = pointer(centre)
  ncnt = length(centre)

  assert(ncnt == dimensions(mesh))

  icxx"""
    // Copy the center location to a Toast pointer
    Point *cnt = NULL;

    if($ncnt == 2)
      cnt = new Point2D($(pcnt)[0], $(pcnt)[1]);

    if($ncnt == 3)
      cnt = new Point3D($(pcnt)[0], $(pcnt)[1], $(pcnt)[2]);

    RVector qm($nvtx, $prhs, SHALLOW_COPY);
    qm = QVec_Point(*$(mesh.ptr), *cnt, SRCMODE_ISOTROPIC);

    delete cnt;
  """

  return rhs

end
