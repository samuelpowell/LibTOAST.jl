# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Assembly" begin

  # Check we can create system matrices and convert to Julia
  S1 = SystemMatrix(mesh3D)
  assemble!(S1, FF)
  @test norm(sparse(S1),1)-norm(sparse(assemble(mesh3D, FF)),1) < 1e-14

  # Check BNDFF integral
  S1 = SystemMatrix(mesh3D)
  assemble!(S1, BNDFF)
  @test norm(sparse(S1),1)-norm(sparse(assemble(mesh3D, BNDFF)),1) < 1e-14

  # Check FF integral
  @test sum(sparse(assemble(mesh3D, FF))) ≈ fullsize(mesh3D)

  # Check DD integral
  @test sum(sparse(assemble(mesh3D, DD))) < 1e-12

  # Check we can add to system matrices in place
  assemble!(S1, DD)
  S2 = sparse(assemble(mesh3D, FF)) + sparse(assemble(mesh3D, DD))
  @test norm(sparse(S1),1)-norm(sparse(S2),1) < 1e-14

  # Check we can assemble with parameters
  p1 = NodalCoeff(mesh3D,rand(nodecount(mesh3D)))
  p2 = NodalCoeff(mesh3D,rand(nodecount(mesh3D)))
  p3 = NodalCoeff(mesh3D,rand(nodecount(mesh3D)))

  S1 = SystemMatrix(mesh3D)
  assemble!(S1, PFF, p1)
  assemble!(S1, PDD, p2)
  assemble!(S1, BNDPFF, p3)

  S2 = sparse(assemble(mesh3D, PFF, p1)) +
       sparse(assemble(mesh3D, PDD, p2)) +
       sparse(assemble(mesh3D, BNDPFF, p3))

  @test norm(sparse(S1),1)-norm(S2,1) < 1e-14

end
