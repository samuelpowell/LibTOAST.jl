# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

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
type PixelMap <: Raster

  ptr::Cxx.CppPtr
  mesh::Mesh

  function PixelMap{N}(mesh::Mesh, bdim::NTuple{N,Integer}; gscale::Integer=2)
    rasterptr = _rasterast_new(mesh, Pixel, bdim, gscale)
    pixelmap = new(rasterptr, mesh)
    finalizer(pixelmap, _raster_delete)
    return pixelmap
  end

end

# Create and initialise a new basis mapper
function _rasterast_new{N}(mesh::Mesh,
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
  dims = Vector{Int}(ndim)
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
