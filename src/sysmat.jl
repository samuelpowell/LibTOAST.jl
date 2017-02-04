# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: sparse

# Export types
export SystemMatrix

immutable SystemMatrix

  ptr::Cxx.CppPtr
  mesh::Mesh

  function SystemMatrix(mesh)

    nnd = numnodes(mesh)
    rowptr = mesh.rowptr
    colidx = mesh.colidx
    matptr = @cxxnew RCompRowMatrix(nnd, nnd, pointer(rowptr), pointer(colidx))

    sysmat = new(matptr, mesh)

    finalizer(sysmat, _delete_sysmat)

    return sysmat

  end

end

# Delete underlying pointer to a Toast++ system matrix
_sysmat_delete(sysmat::SystemMatrix) = icxx"""delete $(sysmat.ptr);"""

function sparse(sysmat::SystemMatrix)
  nnd = numnodes(mesh)
  rowptr = sysmat.mesh.rowptr
  colidx = sysmat.mesh.colidx
  _CSR0_to_CSC1(nnd, nnd, rowptr, colidx, @cxx F.ptr->ValPtr())
end
