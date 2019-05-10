__precompile__(false)

# LibTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2019 Samuel Powell

module libTOAST

using SparseArrays
using LinearAlgebra
using Libdl
using Compat
using Cxx

# Set verbosity, provide info function
verbose = true
_info(s) = verbose && @info(s)

# Load Toast++ libarary path
include(joinpath(dirname(@__FILE__), "..", "deps", "path.jl"))

function __init__()

  # Add header locations: toast, include, libfe, libmath
  addHeaderDir(_jl_toast_dir; kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "include"); kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "src", "libfe"); kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "src", "libmath"); kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "src", "libstoast"); kind = C_System)

  # Define header options accroding to library build options
  defineMacro("TOAST_THREAD")     # Enable threading
  defineMacro("FDOT")             # Enable flourescence (projection)

  # Include headers: felib (this includes all required headers)
  cxxinclude("mathlib.h")         # Matrices and vectors, tasks and verbosity
  cxxinclude("felib.h")           # Mesh related functions
  cxxinclude("stoastlib.h")
  # cxxinclude("source.h")          # Source and detector profiles

  # Import dynamic libraries: libsuperlu, libmath, libfe
  dlopen(_jl_toast_libsuperlu, Libdl.RTLD_GLOBAL)
  dlopen(_jl_toast_libmath, Libdl.RTLD_GLOBAL)
  dlopen(_jl_toast_libfe, Libdl.RTLD_GLOBAL)
  dlopen(_jl_toast_libstoast, Libdl.RTLD_GLOBAL)

  # Initialise Toast++ thread pool
  @cxx Task_Init(0)

end

include("util.jl")            # Utility functions
include("mesh.jl")            # Mesh types
include("sysmat.jl")          # System matrices
include("raster.jl")       # Raster mapping
include("coefftypes.jl")      # Coefficients
include("assembly.jl")        # System matrix assembly
include("regul.jl")  # Regularisation

end # module
