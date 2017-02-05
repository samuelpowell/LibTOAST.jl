# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Meshes" begin

  # Load 2D .msh file from disk and create from vertices, elements
  m = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));
  @test numnodes(m) == 3511
  @test numelems(m) == 6840
  @test numdims(m) == 2
  @test boundingbox(m) == [-25.0 -25.0; 25.0 25.0]
  @test fullsize(m) == 1963.0959349692218

  # Load 3D .msh file from disk and create from vertices, elements
  m = Mesh(joinpath(meshdir_3D, "cyl3.msh"));
  @test numnodes(m) == 27084
  @test numelems(m) == 141702
  @test numdims(m) == 3
  @test boundingbox(m) == [-25.0 25.0; 25.0 -25.0; -25.0 25.0]
  @test fullsize(m) == 98056.71142861326


end
