# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Import
import Base: convert, size, linearindexing, getindex, setindex!, zero, similar

# Export types
export RasterBases, NodalCoeff, SolutionCoeff, RasterCoeff, IntermediateCoeff

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
    return new(mesh,data)
  end
end

NodalCoeff(mesh::Mesh) = NodalCoeff(mesh, Vector{Float64}(nodecount(mesh)))

Base.linearindexing{T<:NodalCoeff}(::Type{T}) = Base.LinearFast()
Base.similar{Te}(coeff::NodalCoeff, ::Type{Te}, dims::Dims) = NodalCoeff(coeff.mesh)
size(coeff::NodalCoeff) = (length(coeff.data),)
getindex(coeff::NodalCoeff, i::Int) = coeff.data[i]
setindex!(coeff::NodalCoeff, v, i::Int) = (coeff.data[i] = v)



# Raster coefficients represent functions in one of the raster bases
abstract RasterCoeffTypes <: AbstractArray{Float64, 1}
Base.linearindexing{T<:RasterCoeffTypes}(::Type{T}) = Base.LinearFast()
Base.similar{T<:RasterCoeffTypes,Te}(coeff::T, ::Type{Te}, dims::Dims) = T(coeff.rmap)
getindex{T<:RasterCoeffTypes}(coeff::T, i::Int) = coeff.data[i]
setindex!{T<:RasterCoeffTypes}(coeff::T, v, i::Int) = (coeff.data[i] = v)

# Solution coefficients represent a function expressed in the solution basis,
# which does not include raster points which fall outside of the support of the
# associated mesh.
type SolutionCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function SolutionCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == slen(rmap))
    return new(rmap,data)
  end
end

size(coeff::SolutionCoeff) = (length(coeff.data),)

function SolutionCoeff(rmap::RasterMap, ci::NodalCoeff)
  co = SolutionCoeff(rmap)
  map!(co,ci)
  return co
end

SolutionCoeff(rmap::RasterMap) = SolutionCoeff(rmap, Vector{Float64}(slen(rmap)))

# Raster coefficients represent a function expressed in the rasterised basis,
# which is defined over a square or cuboid redion, and may include superfluous
# elements which are outside of the support of the mesh.
type RasterCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function RasterCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == blen(rmap))
    return new(rmap,data)
  end
end

size(coeff::RasterCoeff) = (length(coeff.data),)

function RasterCoeff(rmap::RasterMap, ci::NodalCoeff)
  co = RasterCoeff(rmap)
  map!(co,ci)
  return co
end

RasterCoeff(rmap::RasterMap) = RasterCoeff(rmap, Vector{Float64}(blen(rmap)))

# Intermediate coefficients represent a function expressed in a higher resolution
# version of the raster basis, and may include superfluous elements which are
# outside of the support of the mesh.
type IntermediateCoeff <: RasterCoeffTypes
  rmap::RasterMap
  data::Vector{Float64}
  function IntermediateCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == glen(rmap))
    return new(rmap,data)
  end
end

size(coeff::IntermediateCoeff) = (length(coeff.data),)

function IntermediateCoeff(rmap::RasterMap, ci::NodalCoeff)
  co = IntermediateCoeff(rmap)
  map!(co,ci)
  return co
end

IntermediateCoeff(rmap::RasterMap) = IntermediateCoeff(rmap, Vector{Float64}(glen(rmap)))

#
# Coefficient mapping (and construction)
#

# Convert everything to nodal coefficients
function convert{T<:RasterCoeffTypes}(::Type{NodalCoeff}, ci::T)
  co = NodalCoeff(ci.rmap.mesh)
  map!(co, ci)
  return co
end

map!(co::NodalCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rmap, _sn)
map!(co::NodalCoeff, ci::RasterCoeff) = _map!(co.data, ci.data, ci.rmap, _bn)
map!(co::NodalCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rmap, _gn)

# Convert everything to raster coefficients
function convert{T<:RasterCoeffTypes}(::Type{RasterCoeff}, ci::T)
  co = RasterCoeff(ci.rmap)
  map!(co,ci)
  return co
end


map!(co::RasterCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rmap, _nb)
map!(co::RasterCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rmap, _sb)
map!(co::RasterCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rmap, _gb)

# Convert everything to solution coefficients
function convert{T<:RasterCoeffTypes}(::Type{SolutionCoeff}, ci::T)
  co = SolutionCoeff(ci.rmap)
  map!(co,ci)
  return co
end

map!(co::SolutionCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rmap, _ns)
map!(co::SolutionCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, ci.rmap, _bs)
map!(co::SolutionCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rmap, _gs)

# Convert everything to intermediate coefficients
function convert{T<:RasterCoeffTypes}(::Type{IntermediateCoeff}, ci::T)
  co = IntermediateCoeff(ci.rmap)
  map!(co,ci)
  return co
end

map!(co::IntermediateCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rmap, _ng)
map!(co::IntermediateCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, ci.rmap, _bg)
map!(co::IntermediateCoeff, ci::SolutionCoeff) = _map!(co.data, ci.data, ci.rmap, _sg)

# Low level mapping function
function _map!(ovec::Vector{Float64},
               ivec::Vector{Float64},
               rmap::RasterMap,
               mode::MapTransforms)

  modeint = Int(mode)
  ovecptr = pointer(ovec)
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
              $(rmapptr)->Map_MeshToBasis(iprm, oprm);
              break;
          case 1:
              $(rmapptr)->Map_MeshToGrid(iprm, oprm);
              break;
          case 2:
              $(rmapptr)->Map_MeshToSol(iprm, oprm);
              break;

          case 3:
              $(rmapptr)->Map_BasisToMesh(iprm, oprm);
              break;
          case 4:
              $(rmapptr)->Map_BasisToGrid(iprm, oprm);
              break;
          case 5:
              $(rmapptr)->Map_BasisToSol(iprm, oprm);
              break;

          case 6:
              $(rmapptr)->Map_SolToMesh(iprm, oprm);
              break;
          case 7:
              $(rmapptr)->Map_SolToBasis(iprm, oprm);
              break;
          case 8:
              $(rmapptr)->Map_SolToGrid(iprm, oprm);
              break;

          case 9:
              $(rmapptr)->Map_GridToMesh(iprm, oprm);
              break;
          case 10:
              $(rmapptr)->Map_GridToSol(iprm, oprm);
              break;
          case 11:
              $(rmapptr)->Map_GridToBasis(iprm, oprm);
              break;
      }
  """

end
