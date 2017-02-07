# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

@testset "Meshes" begin

  # Load 2D .msh file from disk and create from vertices, elements
  mesh1 = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));
  vtx, ele, elt = data(mesh1)
  mesh2 = Mesh(vtx,ele,elt)

  for m in [mesh1 mesh2]
    @test numnodes(m) == 3511
    @test numelems(m) == 6840
    @test numdims(m) == 2
    @test boundingbox(m) == [-25.0 -25.0; 25.0 25.0]
    @test fullsize(m) == 1963.0959349692218
    @test maxnodes(m) == 3
    vtx, ele = surface(m)
    @test vtx[10,:] == [19.7003,15.3915]
    @test ele[1:5] == [90,135,89,164,88]
  end


  # Load 3D .msh file from disk and create from vertices, elements
  mesh1 = Mesh(joinpath(meshdir_3D, "cyl3.msh"));
  vtx, ele, elt = data(mesh1)
  mesh2 = Mesh(vtx,ele,elt)

  for m in [mesh1 mesh2]
    @test numnodes(m) == 27084
    @test numelems(m) == 141702
    @test numdims(m) == 3
    @test boundingbox(m) == [-25.0 25.0; 25.0 -25.0; -25.0 25.0]
    @test fullsize(m) == 98056.71142861326
    @test maxnodes(m) == 4
    vtx, ele = surface(m)
    @test vtx[10,:] == [-1.29504,6.49008,25.0]
    @test ele[1:5] == [ 2791,1429,3014,2552,124]
  end


end