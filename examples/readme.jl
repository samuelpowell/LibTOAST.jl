# LibTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2019 Samuel Powell

# This example contains a listing of the source code given in the first example
# of the README.md package documentation

# Load the TOAST module
using TOAST

# Load a mesh from disk.
meshdir = joinpath(TOAST._jl_toast_dir, "test", "2D", "meshes")
mesh = Mesh(joinpath(meshdir, "circle25_32.msh"));

# Construct nodal coefficient vectors representing the parameters
μ = fill(NodalCoeff, mesh, 0.01);
κ = fill(NodalCoeff, mesh, 0.3);
ζ = fill(NodalCoeff, mesh, 0.5);

# Build the system matrix in-place
S = SystemMatrix(mesh)
assemble!(S, PDD, κ);
assemble!(S, PFF, μ);
assemble!(S, BNDPFF, ζ);

# Construct a special-purpose source term
source = PointSource(mesh, [0.,0.], Isotropic)
Q = assemble(source)

# Convert the system matrix from its CSR form to a Julia CSC sparse matrix
sysmat = sparse(S)

# Solve as you like
ϕ = sparse(S)\Q

# If you have PyPlot installed...
# using PyPlot
# vtx, = data(mesh);
# plot_trisurf(vtx[:,1], vtx[:,2], ϕ)
