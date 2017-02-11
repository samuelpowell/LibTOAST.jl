# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: string, print, show

# Export types
export RasterMap, PixelMap

# Export methods
export nlen, glen, blen, slen, map, gradient

abstract RasterMap

@enum RasterBases Pixel CubicPixel GaussBlob BesselBlob HanningBlob RampBlob SplineBlob

"""
    PixelMap(mesh, dims; gscale=2)

Map from nodal coefficients to a rasterised grid of pixels of specified
dimension using linear interpolation via an intermediate grid containing an
integer multiple, gscale, of pixels.
"""
type PixelMap <: RasterMap

  ptr::Cxx.CppPtr
  mesh::Mesh

  function PixelMap{N}(mesh::Mesh, bdim::NTuple{N,Integer}; gscale::Integer=2)
    rasterptr = _rastermap_new(mesh, Pixel, bdim, gscale)
    pixelmap = new(rasterptr, mesh)
    finalizer(pixelmap, _raster_delete)
    return pixelmap
  end

end

# Create and initialise a new basis mapper
function _rastermap_new{N}(mesh::Mesh,
                           basis::RasterBases,
                           bdim::NTuple{N,Integer},
                           gscale::Integer;
                           blobarg::Float64 = 1.0,
                           blobrad::Float64 = 1.0)

  bdimarr = [Cint(i) for i in bdim]
  gdimarr = gscale*bdimarr
  dim = dimensions(mesh)
  assert(length(bdimarr)==dim)

  meshptr = mesh.ptr
  bdimptr = pointer(bdimarr)
  gdimptr = pointer(gdimarr)

  basiscode = Int(basis)

  ptr = icxx"""

    RDenseMatrix *bb = 0;

    IVector bdim($dim);
    for (int i = 0; i < $dim; i++)
        bdim[i] = $(bdimptr)[i];

    IVector gdim($(dim));
    for (int i=0; i < $(dim); i++)
        gdim[i] = $(gdimptr)[i];

    Raster *raster = 0;

    switch($basiscode) {
      case 0:
        raster = new Raster_Pixel (bdim, gdim, $meshptr, bb);
        break;
      case 1:
        raster = new Raster_CubicPixel (bdim, gdim, $meshptr, bb);
        break;
      case 2:
        raster = new Raster_GaussBlob (bdim, gdim, $meshptr, $blobarg, $blobrad, bb);
        break;
      case 3:
        raster = new Raster_BesselBlob (bdim, gdim, $meshptr, $blobarg, $blobrad, bb);
        break;
      case 4:
        raster = new Raster_HanningBlob (bdim, gdim, $meshptr, $blobrad, bb);
        break;
      case 5:
        raster = new Raster_RampBlob (bdim, gdim, $meshptr, $blobrad, bb);
        break;
      case 6:
        raster = new Raster_SplineBlob (bdim, gdim, $meshptr, $blobrad, bb);
        break;
    }

    return raster;
  """

  assert(ptr != Ptr{Void}(0x0000000000000000))

  return ptr;

end

# Delete a raster map
_raster_delete(RasterMap::RasterMap) = finalize(RasterMap.ptr)

"""
    nlen(rastermap)

Return the number of nodal coefficients in the associated mesh basis.
"""
nlen(RasterMap::RasterMap) = nodecount(RasterMap.mesh)

"""
    slen(rastermap)

Return the number of coefficients in the raster solution basis, which does
not include raster points outside of the support of the underlying mesh.
"""
slen(RasterMap::RasterMap) = @cxx RasterMap.ptr->SLen()

"""
    blen(rastermap)

Return the number of coefficients in the raster.
"""
blen(RasterMap::RasterMap) = @cxx RasterMap.ptr->BLen()

"""
    glen(rastermap)

Return the number of coefficients in the raster intermediate basis, which is
used to map between bases.
"""
glen(RasterMap::RasterMap) = @cxx RasterMap.ptr->GLen()
