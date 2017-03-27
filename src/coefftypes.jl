# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Import
import Base: convert, size, linearindexing, getindex, setindex!, zero, similar,
  zero, one, fill, \, gradient

# Export types
export RasterBases, NodalCoeff, SolutionCoeff, RasterCoeff, IntermediateCoeff

export gradient

@enum MapTransforms _nb _ng _ns _bn _bg _bs _sn _sb _sg _gn _gs _gb

typealias Reltypes Union{Float64, Complex{Float64}}

#
# Coefficent types
#

# Nodal coefficients represent functions in the mesh basis
"""
    NodalCoeff(mesh, coeff)

Return a nodal coefficient array for a basis defined on the `mesh` initialised
with the `coeff` vector of length equal to the number of nodes in the basis.
The element type of `coeff` may be `Float64`, or `Complex{Float64}`.
"""
type NodalCoeff{T} <: AbstractArray{T, 1}
  mesh::Mesh
  data::Vector{T}
  function (::Type{NodalCoeff}){T<:Reltypes}(mesh::Mesh, data::Vector{T})
    assert(length(data) == nodecount(mesh))
    return new{T}(mesh,data)
  end
end

"""
    NodalCoeff(mesh)

Return an uninitialised nodal coefficient array for a basis defined on the `mesh`.
"""
NodalCoeff(mesh::Mesh) = NodalCoeff(mesh, Float64)
NodalCoeff{T}(mesh::Mesh, ::Type{T}) = NodalCoeff(mesh, Vector{T}(nodecount(mesh)))

Base.linearindexing{T<:NodalCoeff}(::Type{T}) = Base.LinearFast()
Base.similar{T,Te}(coeff::NodalCoeff{T}, ::Type{Te}, dims::Dims) = NodalCoeff(coeff.mesh, T)
size(coeff::NodalCoeff) = (length(coeff.data),)
getindex(coeff::NodalCoeff, i::Int) = coeff.data[i]
setindex!(coeff::NodalCoeff, v, i::Int) = (coeff.data[i] = v)

\(A::SparseMatrixCSC, q::NodalCoeff) = A\q.data

"""
    one(NodalCoeff, mesh)

Return the multiplicative identity element for a basis defined on the `mesh`.
"""
one(::Type{NodalCoeff}, mesh::Mesh) = NodalCoeff(mesh, ones(nodecount(mesh)))

"""
    zero(NodalCoeff, mesh)

Return the additive identity element for a basis defined on the `mesh`.
"""
zero(::Type{NodalCoeff}, mesh::Mesh) = NodalCoeff(mesh, zeros(nodecount(mesh)))

"""
    fill(NodalCoeff, mesh, value)

Return a constant function of `value` for a basis defined on the `mesh`.
"""
function fill{T}(::Type{NodalCoeff}, mesh::Mesh, value::T)
  NodalCoeff(mesh, fill(value, nodecount(mesh)))
end

# Raster coefficients represent functions in one of the raster bases
abstract RasterCoeffTypes{T} <: AbstractArray{T, 1}
Base.linearindexing{T<:RasterCoeffTypes}(::Type{T}) = Base.LinearFast()
getindex{T<:RasterCoeffTypes}(coeff::T, i::Int) = coeff.data[i]
setindex!{T<:RasterCoeffTypes}(coeff::T, v, i::Int) = (coeff.data[i] = v)

"""
    one(RasterCoeffTypes, Raster)

Return the multiplicative identity element for a raster basis defined by `raster`.
"""
one{T<:RasterCoeffTypes}(::Type{T}, rast::Raster) =  T(rast, ones(size(T, rast)...))

"""
    zero(RasterCoeffTypes, mesh)

Return the additive identity element for a raster basis defined by `raster`.
"""
zero{T<:RasterCoeffTypes}(::Type{T}, rast::Raster) = T(rast, zeros(size(T, rast)...))

"""
    fill(RasterCoeffTypes, mesh, value)

Return a constant function of `value` for a basis defined on the `mesh`.
"""
function one{T<:RasterCoeffTypes}(::Type{T}, rast::Raster, value::Float64)
  T(rast, fill(size(T, rast)...))
end

# Solution coefficients represent a function expressed in the solution basis,
# which does not include raster points which fall outside of the support of the
# associated mesh.
"""
    SolutionCoeff(raster, coeff)

Return a coefficient array for the solution basis defined by the `raster`,
initialised with the `coeff` vector of length equal to slen(raster). The element
type of `coeff` may be `Float64`, or `Complex{Float64}`.
"""
type SolutionCoeff{T} <: RasterCoeffTypes{T}
  rast::Raster
  data::Vector{T}
  function (::Type{SolutionCoeff}){T<:Reltypes}(rast::Raster, data::Vector{T})
    assert(length(data) == slen(rast))
    return new{T}(rast,data)
  end
end

Base.similar{T,Te}(coeff::SolutionCoeff{T}, ::Type{Te}, dims::Dims) = SolutionCoeff(coeff.rast, T)
size(coeff::SolutionCoeff) = (length(coeff.data),)
size(::Type{SolutionCoeff}, raster::Raster) = (slen(raster),)

"""
    SolutionCoeff(raster)

Return an uninitialised coefficient array for solution basis defined by the `raster`.
"""
function SolutionCoeff{T}(rast::Raster, ci::NodalCoeff{T})
  co = SolutionCoeff(rast, T)
  map!(co,ci)
  return co
end

SolutionCoeff(rast::Raster) = SolutionCoeff(rast, Float64)
SolutionCoeff{T}(rast::Raster, ::Type{T}) = SolutionCoeff(rast, Vector{T}(slen(rast)))

# Raster coefficients represent a function expressed in the rasterised basis,
# which is defined over a square or cuboid redion, and may include superfluous
# elements which are outside of the support of the mesh.
"""
    RasterCoeff(raster, coeff)

Return a coefficient array for the raster basis defined by the `raster`,
initialised with the `coeff` vector of length equal to blen(raster). The element
type of `coeff` may be `Float64`, or `Complex{Float64}`.
"""
type RasterCoeff{T} <: RasterCoeffTypes{T}
  rast::Raster
  data::Vector{Float64}
  function (::Type{RasterCoeff}){T<:Reltypes}(rast::Raster, data::Vector{T})
    assert(length(data) == blen(rast))
    return new{T}(rast,data)
  end
end

Base.similar{T,Te}(coeff::RasterCoeff{T}, ::Type{Te}, dims::Dims) = RasterCoeff(coeff.rast, T)
size(coeff::RasterCoeff) = (length(coeff.data),)
size(::Type{RasterCoeff}, raster::Raster) = (blen(raster),)

"""
    RasterCoeff(raster)

Return an uninitialised coefficient array for raster basis defined by the `raster`.
"""
function RasterCoeff{T}(rast::Raster, ci::NodalCoeff{T})
  co = RasterCoeff(rast, T)
  map!(co,ci)
  return co
end

RasterCoeff(rast::Raster) = RasterCoeff(rast, Float64)
RasterCoeff{T}(rast::Raster, ::Type{T}) = RasterCoeff(rast, Vector{T}(blen(rast)))

# Intermediate coefficients represent a function expressed in a higher resolution
# version of the raster basis, and may include superfluous elements which are
# outside of the support of the mesh.
type IntermediateCoeff{T} <: RasterCoeffTypes{T}
  rast::Raster
  data::Vector{T}
  function (::Type{IntermediateCoeff}){T<:Reltypes}(rast::Raster, data::Vector{T})
    assert(length(data) == glen(rast))
    return new{T}(rast,data)
  end
end

Base.similar{T,Te}(coeff::IntermediateCoeff{T}, ::Type{Te}, dims::Dims) = RasterCoeff(coeff.rast, T)
size(coeff::IntermediateCoeff) = (length(coeff.data),)
size(::Type{IntermediateCoeff}, raster::Raster) = (glen(raster),)

"""
    IntermediateCoeff(raster, coeff)

Return a coefficient array for the intermediate basis defined by the `raster`,
initialised with the `coeff` vector of length equal to glen(raster). The element
type of `coeff` may be `Float64`, or `Complex{Float64}`.
"""
function IntermediateCoeff{T}(rast::Raster, ci::NodalCoeff{T})
  co = IntermediateCoeff(rast, T)
  map!(co,ci)
  return co
end

"""
    IntermediateCoeff(raster)

Return an uninitialised coefficient array for intermediate basis defined by the `raster`.
"""
IntermediateCoeff(rast::Raster) = IntermediateCoeff(rast, Float64)
IntermediateCoeff{T}(rast::Raster, ::Type{T}) = IntermediateCoeff(rast, Vector{T}(glen(rast)))

#
# Coefficient mapping (and construction)
#

# Convert everything to nodal coefficients
function convert{T<:RasterCoeffTypes}(::Type{NodalCoeff}, ci::T)
  co = NodalCoeff(ci.rast.mesh, eltype(ci.data))
  map!(co, ci)
  return co
end

"""
  map!(out::NodalCoeff, in::RasterCoeffTypes)

Map the function defined by the `input` coefficients defined on a raster basis
to a nodal basis defined on the mesh, overwriting `out` in place.
"""
map!(co::NodalCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rast, _sn)
map!(co::NodalCoeff, ci::RasterCoeff) = _map!(co.data, ci.data, ci.rast, _bn)
map!(co::NodalCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gn)

# Convert everything to raster coefficients
function convert{T<:RasterCoeffTypes}(::Type{RasterCoeff}, ci::T)
  co = RasterCoeff(ci.rast, eltype(ci.data))
  map!(co,ci)
  return co
end

map!(co::RasterCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rast, _nb)
map!(co::RasterCoeff, ci::SolutionCoeff)  = _map!(co.data, ci.data, ci.rast, _sb)
map!(co::RasterCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gb)

# Convert everything to solution coefficients
function convert{T<:RasterCoeffTypes}(::Type{SolutionCoeff}, ci::T)
  co = SolutionCoeff(ci.rast, eltype(ci.data))
  map!(co,ci)
  return co
end

map!(co::SolutionCoeff, ci::NodalCoeff) = _map!(co.data, ci.data, co.rast, _ns)
map!(co::SolutionCoeff, ci::RasterCoeff)  = _map!(co.data, ci.data, ci.rast, _bs)
map!(co::SolutionCoeff, ci::IntermediateCoeff) = _map!(co.data, ci.data, ci.rast, _gs)

# Convert everything to intermediate coefficients
function convert{T<:RasterCoeffTypes}(::Type{IntermediateCoeff}, ci::T)
  co = IntermediateCoeff(ci.rast, eltype(ci.data))
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

# Low level mapping function
function _map!(ovec::Vector{Complex{Float64}},
               ivec::Vector{Complex{Float64}},
               rast::Raster,
               mode::MapTransforms)

  modeint = Int(mode)
  rastptr = rast.ptr

  ilen = length(ivec)

  rivec = real(ivec)
  iivec = imag(ivec)
  rivecptr = pointer(rivec)
  iivecptr = pointer(iivec)

  olen = length(ovec)

  rovec = real(ovec)
  iovec = imag(ovec)
  rovecptr = pointer(rovec)
  iovecptr = pointer(iovec)

  icxx"""

      CVector iprm($(ilen));
      CVector oprm($(olen));

      std::complex<double> *val;

      val = iprm.data_buffer();

      for (int i = 0; i < $(ilen); i++)
          val[i] = std::complex<double> ($(rivecptr)[i], $(iivecptr)[i]);

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

      val = oprm.data_buffer();

      for (int i = 0; i < $(olen); i++)
      {
        $(rovecptr)[i] = real(val[i]);
        $(iovecptr)[i] = imag(val[i]);
      }
  """

  ovec .= rovec .+ im.*iovec

  return ovec

end


"""
    gradient(coeff)

Compute the spatial gradient of `coeff` which must be expressed in an
`IntermediateCoeff` basis.
"""
function gradient(coeff::IntermediateCoeff)

  # TODO: Native implementaiton from ADMM.jl
  ndim = dimensions(coeff.rast.mesh)
  len = glen(coeff.rast)
  rptr = coeff.rast.ptr
  coptr = pointer(coeff.data)

  ∇coeff = Array(Float64, len, ndim)
  gcoptr = pointer(∇coeff)

  icxx"""
    const IVector &gdim = $(rptr)->GDim();
    const RVector &gsize = $(rptr)->GSize();

    RVector img($(len), $(coptr), SHALLOW_COPY);
    RVector *imgrad = new RVector[$(ndim)];
    ImageGradient(gdim, gsize, img, imgrad, $(rptr)->Elref());

    for(int i = 0; i < $(len); i++)
      $(gcoptr)[i] = imgrad[0][i];

    for(int i=0; i < $(len); i++)
      $(gcoptr)[i+$(len)] = imgrad[1][i];

    if($ndim > 2)
    {
      for(int i=0; i < $(len); i++)
        $(gcoptr)[i+(2*$(len))] = imgrad[2][i];
    }

    delete []imgrad;
  """

  return ∇coeff

end

function gradient(coeff::IntermediateCoeff{Complex{Float64}})

  # TODO: Native implementaiton from ADMM.jl
  ndim = dimensions(coeff.rast.mesh)
  len = glen(coeff.rast)
  rptr = coeff.rast.ptr

  rcoeff = real(coeff.data)
  icoeff = imag(coeff.data)
  rcoptr = pointer(rcoeff)
  icoptr = pointer(icoeff)

  ∇coeff = Array(Float64, len*2, ndim)
  gcoptr = pointer(∇coeff)

  icxx"""
    const IVector &gdim = $(rptr)->GDim();
    const RVector &gsize = $(rptr)->GSize();

    CVector img($(len));

    std::complex<double> *val;

    val = img.data_buffer();

    for (int i = 0; i < $(len); i++)
        val[i] = std::complex<double> ($(rcoptr)[i], $(icoptr)[i]);

    CVector *imgrad = new CVector[$(ndim)];
    ImageGradient(gdim, gsize, img, imgrad, $(rptr)->Elref());

    for(int i = 0; i < $(len); i++)
    {
      $(gcoptr)[2*i]   = real(imgrad[0][i]);
      $(gcoptr)[2*i+1] = imag(imgrad[0][i]);
    }

    for(int i=0; i < $(len); i++)
    {
      $(gcoptr)[2*i+(2*$(len))]   = real(imgrad[1][i]);
      $(gcoptr)[2*i+(2*$(len))+1] = imag(imgrad[1][i]);
    }

    if($ndim > 2)
    {
      for(int i=0; i < $(len); i++)
      {
        $(gcoptr)[2*i+(4*$(len))]   = real(imgrad[2][i]);
        $(gcoptr)[2*i+(4*$(len))+1] = real(imgrad[2][i]);
      }
    }

    delete []imgrad;
  """

  ∇coeffout = ∇coeff[1:2:end, :] + im*∇coeff[2:2:end, :]
  return ∇coeffout

end
