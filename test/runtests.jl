using libTOAST
using LinearAlgebra
using Test

# Get some test data from the Toast++ source distribution
meshdir_2D = joinpath(libTOAST._jl_toast_dir, "test", "2D", "meshes")
meshdir_3D = joinpath(libTOAST._jl_toast_dir, "test", "3D", "meshes")


@testset "libTOAST.jl Core" begin

  include("mesh.jl")

  global mesh2D = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));
  global mesh3D = Mesh(joinpath(meshdir_2D, "circle25_32.msh"));

  include("sysmat.jl")
  include("assembly.jl")
  include("raster.jl")
  #include("coefftypes.jl")
  include("regul.jl")

end
