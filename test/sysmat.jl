# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "SystemMatrix" begin

  m = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));

  # Construct an empty system matrix, make sure its empty
  S = SystemMatrix(m)
  @test TOAST.ncols(S) == 3511
  @test TOAST.nrows(S) == 3511

  nvals = length(S.mesh.colidx)
  @test nvals == 24211

  sysmatvals = unsafe_wrap(Array, TOAST.valptr(S), nvals )
  @test all(sysmatvals .== 0)

  # Do some matrix assembly and make sure it comes out the same
  assemble!(S, FF)

  jlsmatvals = sparse(S).nzval
  @test sysmatvals == jlsmatvals

end
