# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Raster" begin

  @testset "PixelMap" begin

    nx, ny = 20, 20
    pixelmap = PixelMap(mesh2D, (nx,ny))
    @test nlen(pixelmap) == nodecount(mesh2D)
    @test blen(pixelmap) == nx*ny
    @test glen(pixelmap) == nx*ny*2*2
    @test slen(pixelmap) == 352

    gscale = 4
    pixelmap = PixelMap(mesh2D, (nx,ny), gscale=gscale)
    @test nlen(pixelmap) == nodecount(mesh2D)
    @test blen(pixelmap) == nx*ny
    @test glen(pixelmap) == nx*ny*gscale*gscale
    @test slen(pixelmap) == 352

  end

end
