# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Export types
export RasterBases NodalCoeff SolutionCoeff RasterCoeff IntermediateCoeff

@enum MapTransforms _nb _ng _ns _bn _bg _bs _sn _sb _sg _gn _gs _gb

#
# Coefficent types
#

# Nodal coefficients represent functions in the mesh basis
type NodalCoeff <: AbstractArray{Float64, 1}
  mesh::Mesh
  data::Vector{Float64}
  function NodalCoeff(mesh::Mesh, data::Vector{Float64})
    assert(length(data) == nodecount(mesh))
  end
end

NodalCoeff(mesh::Mesh) = NodalCoeff(mesh, Vector{Float64}(nodecount(mesh)))

Base.size(coeff::NodalCoeff) = nodecount(coeff.mesh)
Base.linearindexing{T<:NodalCoeff}(::Type{T}) = Base.LinearFast()
Base.getindex(coeff::NodalCoeff, i::Int) = coeff.data[i]

# Raster coefficients represent functions in one of the raster bases
abstract RasterCoeffTypes <: AbstractArray{Float64, 1}
Base.linearindexing{T<:RasterCoeff}(::Type{T}) = Base.LinearFast()
Base.getindex(coeff::RasterCoeff, i::Int) = coeff.data[i]

# Solution coefficients represent a function expressed in the solution basis,
# which does not include raster points which fall outside of the support of the
# associated mesh.
type SolutionCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function SolutionCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == slen(rmap))
  end
end

SolutionCoeff(rmap::RasterMap) = SolutionCoeff(rmap, Vector{Float64}(slen(rmap)))
Base.size(coeff::SolutionCoeff) = slen(coeff.rmap)

# Raster coefficients represent a function expressed in the rasterised basis,
# which is defined over a square or cuboid redion, and may include superfluous
# elements which are outside of the support of the mesh.
type RasterCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function RasterCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == blen(rmap))
  end
end

RasterCoeff(rmap::RasterMap) = RasterCoeff(rmap, Vector{Float64}(blen(rmap)))
Base.size(coeff::RasterCoeff) = blen(coeff.rmap)

# Intermediate coefficients represent a function expressed in a higher resolution
# version of the raster basis, and may include superfluous elements which are
# outside of the support of the mesh.
type IntermediateCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function RasterCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == glen(rmap))
  end
end

IntermediateCoeff(rmap::RasterMap) = IntermediateCoeff(rmap, Vector{Float64}(glen(rmap)))
Base.size(coeff::IntermediateCoeff) = glen(coeff.rmap)

#
# Coefficient mapping (and construction)
#

# Convert everything to nodal coefficients
function convert{T<:RasterCoeffTypes}(NodalCoeff, ci::T)
  co = NodalCoeff(ci.rmap)
  map!(co, ci)
  return co
end

map!(co::NodalCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, co.rmap, _ns)
map!(co::NodalCoeff, ci::RasterCoeff) = _map!((co.data, ci.data, co.rmap, _nb)
map!(co::NodalCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, co.rmap, _ng)

# Convert everything to raster coefficients
function convert{T<:Union{RasterCoeffTypes, NodalCoeff}}(RasterCoeff, ci::T)
  co = RasterCoeffs(ci.rmap)
  map!(co,ci)
  return co
end

map!(co::RasterCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, co.rmap, _sb)
map!(co::RasterCoeff, ci::NodalCoeff) = _map!((co.data, ci.data, co.rmap, _nb)
map!(co::RasterCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, co.rmap, _gb)

# Convert everything to solution coefficients
function convert{T<:Union{RasterCoeffTypes, NodalCoeff}}(SolutionCoeff, ci::T)
  co = SolutionCoeff(ci.rmap)
  map!(co,ci)
  return co
end

map!(co::SolutionCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, co.rmap, _bs)
map!(co::SolutionCoeff, ci::NodalCoeff) = _map!((co.data, ci.data, co.rmap, _ns)
map!(co::SolutionCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, co.rmap, _gs)

# Convert everything to intermediate coefficients
function convert{T<:Union{RasterCoeffTypes, NodalCoeff}}(IntermediateCoeff, ci::T)
  co = IntermediateCoeff(ci.rmap)
  map!(co,ci)
  return co
end

map!(co::IntermediateCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, co.rmap, _bg)
map!(co::IntermediateCoeff, ci::NodalCoeff) = _map!((co.data, ci.data, co.rmap, _ng)
map!(co::IntermediateCoeff, ci::SolutionCoeff) = _map!(co.data, ci.data, co.rmap, _sg)





void raster_map(Raster *raster, double *ivec, int ilen, double *ovec, int olen, int mode)


# Low level mapping function
function _map!(ovec::Vector{Float64},
               ivec::Vector{Float64},
               rmap::RasterMap,
               mode::MapTransforms)

  modeint = Int(mode)
  ovevptr = pointer(ovec)
  ivecptr = pointer(ivec)
  rmapptr = rmap.ptr

  ilen = length(ivec)
  olen = length(ovec)

  icxx"""

      RVector iprm($(ilen), $(ivecptr), SHALLOW_COPY);
      RVector oprm($(olen), $(ovecptr), SHALLOW_COPY);

      switch($(modeint))
      {
          case 0:
              raster->Map_MeshToBasis(iprm, oprm);
              break;
          case 1:
              raster->Map_MeshToGrid(iprm, oprm);
              break;
          case 2:
              raster->Map_MeshToSol(iprm, oprm);
              break;

          case 3:
              raster->Map_BasisToMesh(iprm, oprm);
              break;
          case 4:
              raster->Map_BasisToGrid(iprm, oprm);
              break;
          case 5:
              raster->Map_BasisToSol(iprm, oprm);
              break;

          case 6:
              raster->Map_SolToMesh(iprm, oprm);
              break;
          case 7:
              raster->Map_SolToBasis(iprm, oprm);
              break;
          case 8:
              raster->Map_SolToGrid(iprm, oprm);
              break;

          case 9:
              raster->Map_GridToMesh(iprm, oprm);
              break;
          case 10:
              raster->Map_GridToSol(iprm, oprm);
              break;
          case 11:
              raster->Map_GridToBasis(iprm, oprm);
              break;
      }
  """

end
