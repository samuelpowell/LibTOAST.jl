# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Assembly" begin

  m = Mesh(joinpath(meshdir_3D, "cyl3.msh"))

  # Check we can create system matrices and convert to Julia
  S1 = SystemMatrix(m)
  assemble!(S1, FF)
  @test norm(sparse(S1),1)-norm(sparse(assemble(m, FF)),1) < 1e-14

  # Check BNDFF integral
  S1 = SystemMatrix(m)
  assemble!(S1, BNDFF)
  @test norm(sparse(S1),1)-norm(sparse(assemble(m, BNDFF)),1) < 1e-14

  # Check FF integral
  @test_approx_eq sum(sparse(assemble(m, FF))) fullsize(m)

  # Check DD integral
  @test_skip sum(sparse(assemble(m, DD))) < 1e-12

  # Check we can add to system matrices in place
  assemble!(S1, DD)
  S2 = sparse(assemble(m, FF)) + sparse(assemble(m, DD))
  @test norm(sparse(S1),1)-norm(sparse(S2),1) < 1e-14

  # Check we can assemble with parameters
  p1 = rand(nodecount(m))
  p2 = rand(nodecount(m))
  p3 = rand(nodecount(m))

  S1 = SystemMatrix(m)
  assemble!(S1, PFF, p1)
  assemble!(S1, PDD, p2)
  assemble!(S1, BNDPFF, p3)

  S2 = sparse(assemble(m, PFF, p1)) +
       sparse(assemble(m, PDD, p2)) +
       sparse(assemble(m, BNDPFF, p3))

  @test norm(sparse(S1),1)-norm(S2,1) < 1e-14



end
