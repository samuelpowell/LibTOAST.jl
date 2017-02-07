# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports
import Base: string, print, show

# Export types
export PixelBasis

# Export methods
export map, gradient

# Type definitions
abstract Raster

type PixelBasis <: Raster

  ptr::Cxx.CppPtr
  mesh::Mesh
  #
  function PixelBasis(mesh::Mesh, bdim::Vector{Cint})
  #
  #   gdim = 2*bdim
  #   bb = boundingbobx(mesh)
  #
     ptr = @cxxnew Raster_Pixel
  #   icxx"""$(ptr) = Raster_Pixel ($(pointer(bdim)), $(pointer(gdim)), mesh.ptr, $(pointer(bb)));"""
     raster = new(ptr, mesh)
  #   finalizer(raster, _raster_delete)
  #
  #   return raster
  #
  end

end

_raster_delete(raster::Raster) = finalize(raster.ptr);






#
#
# // raster_new: create a new raster
# Raster * raster_new(Mesh *mesh, int *rdim, int dim, int rtyp, int gscale)
# {
#     int i;
#     RDenseMatrix *bb = 0;
#     double blobrad = 1.0;
#     double blobarg = 1.0;
#     double maptol = 1e-10;
#     double dgscale = 0.1;
#
#     // Copy basis dimensions, set default gdim accordingly
#     IVector bdim(dim);
#     for (i = 0; i < dim; i++)
#         bdim[i] = rdim[i];
#
#     IVector gdim(dim);
#     for (i=0; i <dim; i++)
#         gdim[i] = gscale*bdim[i];
#
#     // Instantiate raster
#     Raster *raster = 0;
#     switch (rtyp) {
#         case 0:
#             raster = new Raster_Pixel (bdim, gdim, mesh, bb);
#             break;
#         case 1:
#             raster = new Raster_CubicPixel (bdim, gdim, mesh, bb);
#             break;
#         case 2:
#             raster = new Raster_GaussBlob (bdim, gdim, mesh, blobarg, blobrad, bb);
#             break;
#         case 3:
#             raster = new Raster_BesselBlob (bdim, gdim, mesh, blobarg, blobrad, bb);
#             break;
#         case 4:
#             raster = new Raster_HanningBlob (bdim, gdim, mesh, blobrad, bb);
#             break;
#         case 5:
#             raster = new Raster_RampBlob (bdim, gdim, mesh, blobrad, bb);
#             break;
#         case 6:
#             raster = new Raster_SplineBlob (bdim, gdim, mesh, blobrad, bb);
#             break;
#         case 7:
#             raster = Raster2::Create<Raster_Pixel2> (bdim, bdim, mesh, bb, maptol);
#             break;
#         case 8:
#             raster = Raster2::Create<Raster_CPixel> (bdim, bdim, mesh, bb, maptol);
#             break;
#         case 9:
#             raster = Raster_Blob2::Create<Raster_Blob2_RB> (bdim, bdim, mesh, blobrad, blobarg, dgscale, bb, maptol);
#             break;
#         case 10:
#             raster = Raster_Blob2::Create<Raster_Blob2_BB> (bdim, bdim, mesh, blobrad, blobarg, dgscale, bb, maptol);
#             break;
#         case 11:
#             raster = Raster_Blob2::Create<Raster_Blob2_SB> (bdim, bdim, mesh, blobrad, blobarg, dgscale, bb, maptol);
#             break;
#         case 12:
#             raster = Raster_Blob2::Create<Raster_Blob2_HB> (bdim, bdim, mesh, blobrad, blobarg, dgscale, bb, maptol);
#             break;
#         case 13:
#             raster = Raster_Blob2::Create<Raster_Blob2_GB> (bdim, bdim, mesh, blobrad, blobarg, dgscale, bb, maptol);
#             break;
#         case 14:
#             raster = Raster2::Create<Raster_CPixel_Tree> (bdim, bdim, mesh, bb, maptol);
#             break;
#     }
#
#     if (bb) delete bb;
#     return raster;
# }
