# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Regularisation" begin

  mesh = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));
  nx, ny = 200, 200
  pixelmap = PixelMap(mesh, (nx,ny))
  vtx,  = data(mesh)

  # Create some reasonably smooth function
  nprm = sin(5*vtx[:,1]./15) + cos(4*vtx[:,2]./15)

  nc = NodalCoeff(mesh, nprm)
  bc = RasterCoeff(pixelmap, nc)

  @testset "TK0" begin

    reg = RegulTK0(bc)

    r0 = RasterCoeff(pixelmap, zeros(blen(pixelmap)))

    @test value(reg, bc) == 0
    @test_approx_eq value(reg, r0) norm(bc.data)^2

    @test norm(TOAST.gradient(reg,bc)) == 0
    @test_approx_eq norm(TOAST.gradient(reg, r0)) 2*norm(bc)

  end

end
