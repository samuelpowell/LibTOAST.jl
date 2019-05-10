# LibTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2019 Samuel Powell

# Imports
import Base: string, print, show

# Export types
export Raster, PixelMap

# Export methods
export nlen, glen, blen, slen, map, bsize, bdim

"""
    Raster

A rasterisation of the bounding box of a geometry defined by the mesh, with
methods capable of transforming functions between nodal, and rasterised bases.

See PixelMap.
"""
@compat abstract type Raster end

@enum RasterBases Pixel CubicPixel GaussBlob BesselBlob HanningBlob RampBlob SplineBlob

"""
    PixelMap(mesh, dims; gscale=2)

Map from nodal coefficients to a rasterised grid of pixels of specified
dimension using linear interpolation via an intermediate grid containing an
integer multiple, gscale, of pixels.

# Arguments
* `mesh::Mesh`: a TOAST mesh
* `dims::NTuple{ndim,Integer}`: the dimensions of the raster basis
* `gscale=2::Int`: the scale factor used for the intermediate mapping raster
"""
mutable struct PixelMap <: Raster

  ptr::Cxx.CxxCore.CppPtr
  mesh::Mesh

  function PixelMap(mesh::Mesh, bdim::NTuple{N,Integer}; gscale::Integer=2) where {N}
    rasterptr = _rasterast_new(mesh, Pixel, bdim, gscale)
    pixelmap = new(rasterptr, mesh)
    finalizer(_raster_delete, pixelmap)
    return pixelmap
  end

end

# Create and initialise a new basis mapper
function _rasterast_new(mesh::Mesh,
                        basis::RasterBases,
                        dims::NTuple{N,Integer},
                        gscale::Integer;
                        blobarg::Float64 = 1.0,
                        blobrad::Float64 = 1.0) where {N}

  bdimarr = [Cint(i) for i in dims]
  gdimarr = [Cint(i*gscale) for i in dims]
  basiscode = Cint(basis)
  dim = Cint(dimensions(mesh))

  @assert length(bdimarr)==dim

  meshptr = mesh.ptr

  bd = @cxxnew IVector(Cint(2), pointer(bdimarr), icxx"""return DEEP_COPY;""")
  gd = @cxxnew IVector(Cint(2), pointer(gdimarr), icxx"""return DEEP_COPY;""")
  
  ptr = icxx"""

    RDenseMatrix *bb = 0;
    Raster *raster = 0;

    switch($basiscode) {
      case 0:
        raster = new Raster_Pixel (*$(bd), *$(gd), $meshptr, bb);
        break;
      case 1:
        raster = new Raster_CubicPixel (*$(bd), *$(gd), $meshptr, bb);
        break;
      case 2:
        raster = new Raster_GaussBlob (*$(bd),*$(gd), $meshptr, $blobarg, $blobrad, bb);
        break;
      case 3:
        raster = new Raster_BesselBlob (*$(bd), *$(gd), $meshptr, $blobarg, $blobrad, bb);
        break;
      case 4:
        raster = new Raster_HanningBlob (*$(bd), *$(gd), $meshptr, $blobrad, bb);
        break;
      case 5:
        raster = new Raster_RampBlob (*$(bd), *$(gd), $meshptr, $blobrad, bb);
        break;
      case 6:
        raster = new Raster_SplineBlob (*$(bd), *$(gd), $meshptr, $blobrad, bb);
        break;
    }

    return raster;
  """

  @assert ptr != C_NULL

  return ptr;

end

# Delete a raster map
_raster_delete(Raster::Raster) = finalize(Raster.ptr)

"""
    nlen(raster)

Return the number of nodal coefficients in the mesh basis associated with the
raster.
"""
nlen(Raster::Raster) = nodecount(Raster.mesh)

"""
    slen(raster)

Return the number of coefficients in the raster solution basis, which does
not include raster points outside of the support of the underlying mesh.
"""
slen(Raster::Raster) = @cxx Raster.ptr->SLen()

"""
    blen(raster)

Return the number of coefficients in the raster.
"""
blen(Raster::Raster) = @cxx Raster.ptr->BLen()

"""
    glen(raster)

Return the number of coefficients in the raster intermediate basis, which is
used to map between bases.
"""
glen(Raster::Raster) = @cxx Raster.ptr->GLen()

"""
    bdim(raster)

Return the number of elements in each dimension of the raster basis.
"""
function bdim(raster::Raster)

  ndim = dimensions(raster.mesh)
  dims = Vector{Int}(undef, ndim)
  dimp = pointer(dims)

  icxx"""
    for(int i=0; i<$(ndim); i++)
      $(dimp)[i] = $(raster.ptr)->BDim()[i];
  """

  return dims

end

"""
    bsize(raster)

The `size` (area/volume) of each element of the raster basis.
"""
function bsize(raster::Raster)
  bb = boundingbox(raster.mesh)
  prod((bb[2,:] .- bb[1,:])./ bdim(raster) )
end
