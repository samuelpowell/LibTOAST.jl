# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

module libtoast

using Cxx

# Load Toast++ libarary path
include(joinpath(dirname(@__FILE__), "..", "deps", "path.jl"))

function __init__()

  # Add header locations: toast, include, libfe, libmath
  addHeaderDir(_jl_toast_dir; kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "include"); kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "src", "libfe"); kind = C_System)
  addHeaderDir(joinpath(_jl_toast_dir, "src", "libmath"); kind = C_System)

  # Include headers: felib (this includes all required headers)
  cxxinclude("felib.h");

  # Import dynamic libraries: libsuperlu, libmath, libfe
  Libdl.dlopen(_jl_toast_libsuperlu)
  Libdl.dlopen(_jl_toast_libmath)
  Libdl.dlopen(_jl_toast_libfe)

end

end # module
