using TOAST
using Base.Test

# Get some test data from the Toast++ source distribution
meshdir_2D = joinpath(TOAST._jl_toast_dir, "test", "2D", "meshes")
meshdir_3D = joinpath(TOAST._jl_toast_dir, "test", "3D", "meshes")

include("mesh.jl")
include("sysmat.jl")
include("assembly.jl")
