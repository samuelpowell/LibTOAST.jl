# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# This example contains a listing of the source code given in the first example
# of the README.md package documentation

using TOAST

meshdir = joinpath(TOAST._jl_toast_dir, "test", "2D", "meshes")
mesh = Mesh(joinpath(meshdir, "circle25_32.msh"));

μ = fill(NodalCoeff, mesh, 0.01);
κ = fill(NodalCoeff, mesh, 0.3);
ζ = fill(NodalCoeff, mesh, 0.5);

S = SystemMatrix(mesh)
assemble!(S, PDD, κ);
assemble!(S, PFF, μ);
assemble!(S, BNDPFF, ζ);

vtx, = data(mesh)
source = NodalCoeff(mesh, 1./((vtx[:,1].^2 + vtx[:,2].^2) + 0.1 ))

Q = assemble(mesh, PointSource, [0.,0.])

sysmat = sparse(S)

ϕ = sparse(S)\Q
