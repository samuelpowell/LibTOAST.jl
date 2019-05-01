# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: sparse

# Export types
export SystemMatrix

mutable struct SystemMatrix

  ptr::Cxx.CxxCore.CppPtr
  mesh::Mesh

  function SystemMatrix(mesh)

    nnd = nodecount(mesh)
    rowptr = mesh.rowptr
    colidx = mesh.colidx
    matptr = @cxxnew RCompRowMatrix(nnd, nnd, pointer(rowptr), pointer(colidx))

    sysmat = new(matptr, mesh)

    finalizer(sysmat, _sysmat_delete)

    return sysmat

  end

end

ncols(S::SystemMatrix) = @cxx S.ptr->nCols()
nrows(S::SystemMatrix) = @cxx S.ptr->nRows()
valptr(S::SystemMatrix) = @cxx S.ptr->ValPtr()

# Delete underlying pointer to a Toast++ system matrix
# TODO: Check if Cxx.jl does this for me
_sysmat_delete(sysmat::SystemMatrix) = finalize(sysmat.ptr) #icxx"""delete $(sysmat.ptr);"""

function sparse(sysmat::SystemMatrix)
  nnd = nodecount(sysmat.mesh)
  rowptr = sysmat.mesh.rowptr
  colidx = sysmat.mesh.colidx
  _CSR0_to_CSC1(nnd, nnd, rowptr, colidx, valptr(sysmat))
end
