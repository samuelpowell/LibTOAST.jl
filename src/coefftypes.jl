# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Import
import Base: convert, size, linearindexing, getindex, setindex!, zero, similar,
  zero, one

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
one(::Type{NodalCoeff}, mesh::Mesh) = NodalCoeff(mesh, ones(nodecount(mesh)))
zero(::Type{NodalCoeff}, mesh::Mesh) = NodalCoeff(mesh, zeros(nodecount(mesh)))

# Raster coefficients represent functions in one of the raster bases
abstract RasterCoeffTypes <: AbstractArray{Float64, 1}
Base.linearindexing{T<:RasterCoeffTypes}(::Type{T}) = Base.LinearFast()
Base.similar{T<:RasterCoeffTypes,Te}(coeff::T, ::Type{Te}, dims::Dims) = T(coeff.rast)
getindex{T<:RasterCoeffTypes}(coeff::T, i::Int) = coeff.data[i]
setindex!{T<:RasterCoeffTypes}(coeff::T, v, i::Int) = (coeff.data[i] = v)
function one{T<:RasterCoeffTypes}(::Type{T}, rast::Raster)
  coeff = T(rast)
  coeff .= 1
end
function zero{T<:RasterCoeffTypes}(::Type{T}, rast::Raster)
  coeff = T(rast)
  coeff .= 0
end


# Solution coefficients represent a function expressed in the solution basis,
# which does not include raster points which fall outside of the support of the
# associated mesh.
type SolutionCoeff <: RasterCoeffTypes
  rast::Raster
  data::Vector{Float64}
  function SolutionCoeff(rast::Raster, data::Vector{Float64})
    assert(length(data) == slen(rast))
    return new(rast,data)
  end
end

size(coeff::SolutionCoeff) = (length(coeff.data),)

function SolutionCoeff(rast::Raster, ci::NodalCoeff)
  co = SolutionCoeff(rast)
  map!(co,ci)
  return co
end

SolutionCoeff(rast::Raster) = SolutionCoeff(rast, Vector{Float64}(slen(rast)))

# Raster coefficients represent a function expressed in the rasterised basis,
# which is defined over a square or cuboid redion, and may include superfluous
# elements which are outside of the support of the mesh.
type RasterCoeff <: RasterCoeffTypes
  rast::Raster
  data::Vector{Float64}
  function RasterCoeff(rast::Raster, data::Vector{Float64})
    assert(length(data) == blen(rast))
    return new(rast,data)
  end
end

size(coeff::RasterCoeff) = (length(coeff.data),)

function RasterCoeff(rast::Raster, ci::NodalCoeff)
  co = RasterCoeff(rast)
  map!(co,ci)
  return co
end

RasterCoeff(rast::Raster) = RasterCoeff(rast, Vector{Float64}(blen(rast)))

# Intermediate coefficients represent a function expressed in a higher resolution
# version of the raster basis, and may include superfluous elements which are
# outside of the support of the mesh.
type IntermediateCoeff <: RasterCoeffTypes
  rast::Raster
  data::Vector{Float64}
  function IntermediateCoeff(rast::Raster, data::Vector{Float64})
    assert(length(data) == glen(rast))
    return new(rast,data)
  end
end

size(coeff::IntermediateCoeff) = (length(coeff.data),)

function IntermediateCoeff(rast::Raster, ci::NodalCoeff)
  co = IntermediateCoeff(rast)
  map!(co,ci)
  return co
end

IntermediateCoeff(rast::Raster) = IntermediateCoeff(rast, Vector{Float64}(glen(rast)))

#
# Coefficient mapping (and construction)
#

# Convert everything to nodal coefficients
function convert{T<:RasterCoeffTypes}(::Type{NodalCoeff}, ci::T)
  co = NodalCoeff(ci.rast.mesh)
  map!(co, ci)
  return co
end

map!(co::NodalCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rast, _sn)
map!(co::NodalCoeff, ci::RasterCoeff) = _map!(co.data, ci.data, ci.rast, _bn)
map!(co::NodalCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gn)

# Convert everything to raster coefficients
function convert{T<:RasterCoeffTypes}(::Type{RasterCoeff}, ci::T)
  co = RasterCoeff(ci.rast)
  map!(co,ci)
  return co
end


map!(co::RasterCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rast, _nb)
map!(co::RasterCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rast, _sb)
map!(co::RasterCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gb)

# Convert everything to solution coefficients
function convert{T<:RasterCoeffTypes}(::Type{SolutionCoeff}, ci::T)
  co = SolutionCoeff(ci.rast)
  map!(co,ci)
  return co
end

map!(co::SolutionCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rast, _ns)
map!(co::SolutionCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, ci.rast, _bs)
map!(co::SolutionCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gs)

# Convert everything to intermediate coefficients
function convert{T<:RasterCoeffTypes}(::Type{IntermediateCoeff}, ci::T)
  co = IntermediateCoeff(ci.rast)
  map!(co,ci)
  return co
end

map!(co::IntermediateCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rast, _ng)
map!(co::IntermediateCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, ci.rast, _bg)
map!(co::IntermediateCoeff, ci::SolutionCoeff) = _map!(co.data, ci.data, ci.rast, _sg)

# Low level mapping function
function _map!(ovec::Vector{Float64},
               ivec::Vector{Float64},
               rast::Raster,
               mode::MapTransforms)

  modeint = Int(mode)
  ovecptr = pointer(ovec)
  ivecptr = pointer(ivec)
  rastptr = rast.ptr

  ilen = length(ivec)
  olen = length(ovec)

  icxx"""

      RVector iprm($(ilen), $(ivecptr), SHALLOW_COPY);
      RVector oprm($(olen), $(ovecptr), SHALLOW_COPY);

      switch($(modeint))
      {
          case 0:
              $(rastptr)->Map_MeshToBasis(iprm, oprm);
              break;
          case 1:
              $(rastptr)->Map_MeshToGrid(iprm, oprm);
              break;
          case 2:
              $(rastptr)->Map_MeshToSol(iprm, oprm);
              break;

          case 3:
              $(rastptr)->Map_BasisToMesh(iprm, oprm);
              break;
          case 4:
              $(rastptr)->Map_BasisToGrid(iprm, oprm);
              break;
          case 5:
              $(rastptr)->Map_BasisToSol(iprm, oprm);
              break;

          case 6:
              $(rastptr)->Map_SolToMesh(iprm, oprm);
              break;
          case 7:
              $(rastptr)->Map_SolToBasis(iprm, oprm);
              break;
          case 8:
              $(rastptr)->Map_SolToGrid(iprm, oprm);
              break;

          case 9:
              $(rastptr)->Map_GridToMesh(iprm, oprm);
              break;
          case 10:
              $(rastptr)->Map_GridToSol(iprm, oprm);
              break;
          case 11:
              $(rastptr)->Map_GridToBasis(iprm, oprm);
              break;
      }
  """

end
