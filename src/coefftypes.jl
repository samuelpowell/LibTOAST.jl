# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Export types
export RasterBases NodalCoeff SolutionCoeff RasterCoeff IntermediateCoeff


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
abstract RasterCoeff <: AbstractArray{Float64, 1}
Base.linearindexing{T<:RasterCoeff}(::Type{T}) = Base.LinearFast()
Base.getindex(coeff::RasterCoeff, i::Int) = coeff.data[i]

# Solution coefficients represent a function expressed in the solution basis,
# which does not include raster points which fall outside of the support of the
# associated mesh.
type SolutionCoeff <: RasterCoeff
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
type RasterCoeff <: RasterCoeff
  rmap::RasterMap
  data::Vector{Float64}
  function RasterCoeff(rmap::RasterMap, data::Vector{Float64})
    assert(length(data) == blen(rmap))
  end
end

RasterCoeff(rmap::RasterMap) = RasterCoeff(rmap, Vector{Float64}(blen(rmap)))

Base.size(coeff::RasterCoeff) = blen(coeff.rmap)





function convert{T<:RasterCoeff}(NodalCoeff, ci::T)
  co = NodalCoeff(ci.rmap, Vector{Float64}(nlen(ci.rmap)))
  map!(co, ci)
  return co
end

map!(co::NodalCoeff, ci::SolutionCoeff) assert(co.rmap == ci.rmap)
  _map!(

   = _map!()
  rasterptr


end

function map!(co::NodalCoeff, ci::RasterCoeff)
end

# Convert to solution
function convert(SolutionCoeff, coeff::NodalCoeff)

  = map(SolutionCoeff, coeff)

# Construct a nodal coefficient vector from a raster coefficient vector
function NodalCoeff{T<:RasterCoeff}(coeff::T)
  ncoeff = Vector{Float64}(nlen(coeff.raster))
  map!(rmap, coeff, ncoeff)
  return ncoeff
end

# Construct a coefficient vector in the rasterised basis from
function RasterCoeff(coeff::SolutionCoeff)
  rcoeff = RasterCoeff(rmap.mesh, Vector{Float64}(blen(rmap)))
  map!(rmap, coeff, rcoeff)
  return rcoeff
end

RasterCoeff(raster::RasterCoeff) =

function IntermediateCoeff(raster::RasterCoeff,  coeff::SolutionCoeff)
  icoeff = IntermediateCoeff(rmap.mesh, Vector{Float64}(glen(rmap)))
  map!(rmap, coeff, icoeff)
  return icoeff
end

IntermediateCoeff(raster::RasterCoeff) =

SolutionCoeff(raster::RasterCoeff) =

SolutionCoeff(raster::RasterCoeff) =




# Low level mapping function
function _map!(ovec::Vector{Float64},
               ivec::Vector{Float64},
               mesh::Mesh,
               rmap::RasterMap,
               mode::MapModes)







end
