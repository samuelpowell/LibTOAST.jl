# LibTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2019 Samuel Powell

@testset "Regularisation" begin

  fdiskip = 300     # Finite difference applied to every xth element

  nx, ny = 300, 300
  pixelmap = PixelMap(mesh2D, (nx,ny))
  vtx,  = data(mesh2D)

  # Create some reasonably smooth function
  nprm = sin.(5.*vtx[:,1]./15) .+ cos.(4.*vtx[:,2]./15)

  nc = NodalCoeff(mesh2D, nprm)
  rc = SolutionCoeff(pixelmap, nc)
  rcpdx = similar(rc)
  rcmdx = similar(rc)
  rcp1 = copy(rc)
  rcp1 .+= 1
  r0 = SolutionCoeff(pixelmap, zeros(slen(pixelmap)))

  @testset "TK0" begin

    reg = RegulTK0(rc)

    @test val(reg, rc) == 0
    @test val(reg, r0) ≈ norm(rc.data)^2

    @test norm(libTOAST.grad(reg,rc)) == 0
    @test norm(libTOAST.grad(reg, r0)) ≈ 2*norm(rc)

  end

  @testset "TK1" begin

    reg = RegulTK1(rc)

    @test val(reg, rc) == 0.0
    @test ≈(val(reg, rcp1), 0.0, atol=1e-20)

    # Finite difference the grad
    g = grad(reg, rc)

    for i = 1:fdiskip:length(g)
      dx = 0.0001*rc[i]
      rcpdx .= rc; rcpdx[i] += dx
      rcmdx .= rc; rcmdx[i] -= dx
      @test abs((val(reg, rcpdx) - val(reg, rcmdx))/(2*dx) - g[i]) < 1e-5
    end

  end

  @testset "TV" begin

    reg = RegulTV(rc)

    @test val(reg, rc) == 0.0
    @test ≈(val(reg, rcp1), 0.0, atol=1e-20)

    # Finite difference the grad
    g = grad(reg, rc)

    for i = 1:fdiskip:length(g)
      dx = 0.0001*rc[i]
      rcpdx .= rc; rcpdx[i] += dx
      rcmdx .= rc; rcmdx[i] -= dx
      @test abs((val(reg, rcpdx) - val(reg, rcmdx))/(2*dx) - g[i]) < 1e-5
    end

  end

  @testset "PM" begin

    reg = RegulPM(rc)

    @test val(reg, rc) == 0.0
    @test ≈(val(reg, rcp1), 0.0, atol=1e-20)

    # Finite difference the grad
    g = grad(reg, rc)

    for i = 1:fdiskip:length(g)
      dx = 0.0001*rc[i]
      rcpdx .= rc; rcpdx[i] += dx
      rcmdx .= rc; rcmdx[i] -= dx
      @test abs((val(reg, rcpdx) - val(reg, rcmdx))/(2*dx) - g[i]) < 1e-5
    end

  end

end
