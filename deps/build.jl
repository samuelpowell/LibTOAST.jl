# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Directories
toastv = "v2.0.1"
prefix = dirname(@__FILE__)
mkpath(prefix)

# Get binaries
if is_apple()
  osroot = "darwin"
  suffix = ".zip"
end

if is_windows()
  osroot = "windows"
  suffix = ".zip"
end

if is_linux()
  osroot = "linux"
  suffix = ".tar.gz"
end

osroot = (is_apple() ? "darwin" : (is_linux() ? "linux" : "windows"))
dlroot = "https://github.com/toastpp/toastpp/releases/download/" * toastv * "/"
dlfile = "toast_bin_" * osroot * "64" * suffix

download(dlroot * dlfile, joinpath(prefix, dlfile))

# Unzip
if is_linux()
  run(`tar xzvf $(joinpath(prefix, dlfile)) -C $prefix`)
else
  run(`unzip -o $(joinpath(prefix, dlfile)) -d $prefix`)
end

# Remove temporary download
rm(joinpath(prefix, dlfile))
