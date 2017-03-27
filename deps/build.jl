# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Directories
toastn = "2.0.2"
toastv = "v" * toastn
prefix = dirname(@__FILE__)
mkpath(prefix)

# Get source and binaries
if is_apple()
  osroot = "darwin"
  oslibn = ".dylib"
end

if is_windows()
  osroot = "windows"
  oslibn = ".dll"
end

if is_linux()
  osroot = "linux"
  oslibn = ".so"
end

suffix = ".zip"

srcroot = "https://github.com/toastpp/toastpp/archive/"
srcfile = toastv * suffix
srcdir = "toastpp-" * toastn

info("Downloading Toast++ source")
download(srcroot * srcfile, joinpath(prefix, srcfile))

osroot = (is_apple() ? "darwin" : (is_linux() ? "linux" : "windows"))
dlroot = "https://github.com/toastpp/toastpp/releases/download/" * toastv * "/"
dlfile = "toast_bin_" * osroot * "64" * suffix

info("Downloading Toast++ library binaries")
download(dlroot * dlfile, joinpath(prefix, dlfile))

# Unzip
info("Uncompressing Toast++ source and binaries")
if is_linux()
  run(`unzip -o $(joinpath(prefix, dlfile)) -d $prefix`)
  run(`unzip -o $(joinpath(prefix, srcfile)) -d $prefix`)
else
  run(`unzip -o $(joinpath(prefix, dlfile)) -d $prefix`)
  run(`unzip -o $(joinpath(prefix, srcfile)) -d $prefix`)
end

if is_linux()
  run(`cp -r $(joinpath(prefix, osroot * "64" )) $(joinpath(prefix, "toastpp-" * toastn))`)
  run(`rm -rf $(joinpath(prefix, osroot * "64"))`)
else
  run(`cp -r $(joinpath(prefix, "toast", osroot * "64" )) $(joinpath(prefix, "toastpp-" * toastn))`)
  run(`rm -rf $(joinpath(prefix, "toast"))`)
end

# Remove temporary download
info("Deleteing temporary files")
rm(joinpath(prefix, dlfile))
rm(joinpath(prefix, srcfile))

# Write path
info("Storing Toast++ installation path")
libprefix = joinpath(prefix, "toastpp-" * toastn, osroot * "64", "lib");

s = """
const _jl_toast_dir = "$(joinpath(prefix, "toastpp-" * toastn))"
const _jl_toast_libdir =  "$libprefix"
const _jl_toast_libfe = "$(joinpath(libprefix, "libfe" * oslibn))"
const _jl_toast_libmath = "$(joinpath(libprefix, "libmath" * oslibn))"
const _jl_toast_libstoast = "$(joinpath(libprefix, "libstoast" * oslibn))"
const _jl_toast_libsuperlu = "$(joinpath(libprefix, "libsuperlu" * oslibn))"
"""

f = open(joinpath(prefix, "path.jl"), "w")
write(f, s)
close(f)

# Fin
info("Installation complete")
