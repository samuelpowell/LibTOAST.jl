# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Directories
toastn = "2.0.2"
toastv = "v" * toastn
prefix = dirname(@__FILE__)
mkpath(prefix)

# Get source and binaries
if Sys.isapple()
  osroot = "darwin"
  oslibn = ".dylib"
end

if Sys.iswindows()
  osroot = "windows"
  oslibn = ".dll"
end

if Sys.islinux()
  osroot = "linux"
  oslibn = ".so"
end

suffix = ".zip"

srcroot = "https://github.com/toastpp/toastpp/archive/"
srcfile = toastv * suffix
srcdir = "toastpp-" * toastn

@info "Downloading Toast++ source"
download(srcroot * srcfile, joinpath(prefix, srcfile))

osroot = (Sys.isapple() ? "darwin" : Sys.islinux() ? "linux" : "windows")
dlroot = "https://github.com/toastpp/toastpp/releases/download/" * toastv * "/"
dlfile = "toast_bin_" * osroot * "64" * suffix

@info "Downloading Toast++ library binaries"
download(dlroot * dlfile, joinpath(prefix, dlfile))

# Unzip
@info "Uncompressing Toast++ source and binaries"
if Sys.islinux()
  # On linux, unzip binaries into toast subdirectory to match OS X
  run(`unzip -q -o $(joinpath(prefix, dlfile)) -d $(joinpath(prefix, "toast"))`)
else
  run(`unzip -q -o $(joinpath(prefix, dlfile)) -d $prefix`)
end
run(`unzip -q -o $(joinpath(prefix, srcfile)) -d $prefix`)

# Move the binaries into the source tree
@info "Adding binaries to source tree"
run(`cp -r $(joinpath(prefix, "toast", osroot * "64" )) $(joinpath(prefix, "toastpp-" * toastn))`)
run(`rm -rf $(joinpath(prefix, osroot * "64"))`)

# Remove temporary download
@info "Deleteing temporary files"
rm(joinpath(prefix, dlfile))
rm(joinpath(prefix, srcfile))

# Write path
@info "Storing Toast++ installation path"
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
@info "Installation complete"
