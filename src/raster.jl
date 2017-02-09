# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: string, print, show

# Export types
export PixelRaster

# Export methods
export map, gradient

@enum RasterShapes Pixel CubicPixel GaussBlob BesselBlob HanningBlob RampBlob SplineBlob

# Type definitions
type Raster{S<:RasterShapes}

  ptr::Cxx.CppPtr
  mesh::Mesh

  function Raster{T<:Integer}(shape::RasterShapes, mesh::Mesh, rdim::Vector{T}; gscale=2)

    dim = dimensions(mesh)
    assert(length(rdim)==dim)
    meshptr = mesh.ptr
    rdimptr = pointer(rdim)

    ptr = icxx"""

      RDenseMatrix *bb = 0;

      IVector bdim($(dim));
      for (int i = 0; i < $(dim); i++)
          bdim[i] = $(rdimptr)[i];

      IVector gdim($(dim));
      for (int i=0; i < $(dim); i++)
          gdim[i] = $(gscale)*bdim[i];

      Raster *raster;

      switch($(Int(shape))) {
        case 1:
          raster = new Raster_CubicPixel (bdim, gdim, mesh, bb);
          break;
        case 2:
          raster = new Raster_GaussBlob (bdim, gdim, mesh, blobarg, blobrad, bb);
          break;
        case 3:
          raster = new Raster_BesselBlob (bdim, gdim, mesh, blobarg, blobrad, bb);
          break;
        case 4:
          raster = new Raster_HanningBlob (bdim, gdim, mesh, blobrad, bb);
          break;
        case 5:
          raster = new Raster_RampBlob (bdim, gdim, mesh, blobrad, bb);
          break;
        case 6:
          raster = new Raster_SplineBlob (bdim, gdim, mesh, blobrad, bb);
          break;
      }

      return raster;

    """
    raster = new{shape}(ptr, mesh)
    finalizer(raster, _raster_delete)

    return raster

  end

end

_raster_delete(raster::Raster) = finalize(raster.ptr);

nlen(raster::Raster) = nodecount(raster.mesh)
slen(raster::Raster) = @cxx raster.ptr->SLen()
blen(raster::Raster) = @cxx raster.ptr->Blen()
glen(raster::Raster) = @cxx raster.ptr->Glen()


#
#
# function map(r::Raster, v::Raster, ::rbasis)
#   ovec = vec
#   map!(r::Rater, v::Raster, o::Raster)
# end
#
# function map!(r::Raster, v::Raster, ::rbasis)
#   icxx"""
#     RVector iprm(ilen, ivec, SHALLOW_COPY);
#     RVector oprm(olen, ovec, SHALLOW_COPY);
#     raster->Map_MeshToBasis(iprm, oprm);
#   """
# end
#
