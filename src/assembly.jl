# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Export integral enumerations
export FF, DD, BNDFF, PFF, PDD, BNDPFF, P, PF, BNDPFF

# Methods
export assemble, assemble!

# The following enumerations are defined in mesh.h
@enum BilinearIntegrals FF=0 DD=1 BNDFF=12
@enum BilinearParamIntegrals PFF=2 PDD=3 BNDPFF=4
@enum LinearParamIntegrals P=0 PF=1 BNDPF=2

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
