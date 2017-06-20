# libTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports

# Export types
export Regul, RegulTK0, RegulTK1, RegulTV, RegulPM

# Export methods
export val, grad, hess

"""
    Regul

Abstract supertype of regulariers capable of calculating the value, gradient,
and Hessian of of a particular regularisation functional, with parameters
defined in a SolutionCoeff raster basis.

See RegulTK0, RegulTK1, RegulPM, RegulTV
"""
@compat abstract type Regul end

"""
  RegulTK0(x₀)

Construct a zeroth-order Tikhonov regularisation functional with baseline
parameters x₀ (defined in a SolutionCoeff basis).

``f(x) = ‖x-x₀‖²``

"""
type RegulTK0 <: Regul

  ptr::Cxx.CppPtr
  rast::Raster

  function RegulTK0(x0::SolutionCoeff)

    rast = x0.rast
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      RVector xs($(xlen), 1);
      Regularisation *reg = new Tikhonov0 (1.0, &x0, &xs);
      return reg;
    """

    regul = new(regulptr, rast)
    finalizer(regul, _regul_delete)
    return regul
  end

end

"""
  RegulTK1(x₀)

Construct a first-order Tikhonov regularisation functional with baseline
parameters x₀ (defined in a SolutionCoeff basis).

``f(x) = ‖∇(x-x₀)‖²``
"""
type RegulTK1 <: Regul

  ptr::Cxx.CppPtr
  rast::Raster

  function RegulTK1(x0::SolutionCoeff)

    rast = x0.rast
    rptr = rast.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      void *kapref = 0;
      bool istensor = false;
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new TK1 (1.0, &x0, $(rptr), kapref, istensor);
      return reg;
    """

    regul = new(regulptr, rast)
    finalizer(regul, _regul_delete)
    return regul
  end

end

"""
    RegulTV(x₀, β=1.0)

Construct a soft Total-Variation regularisation functional with parameter
β and baseline parameters x₀ (defined in a SolutionCoeff basis).
"""
type RegulTV <: Regul

  ptr::Cxx.CppPtr
  rast::Raster

  function RegulTV(x0::SolutionCoeff, β::Float64)

    rast = x0.rast
    rptr = rast.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new TV(1.0, $(β), &x0, $(rptr));
      return reg;
    """

    regul = new(regulptr, rast)
    finalizer(regul, _regul_delete)
    return regul

  end

end

function RegulTV(x0::SolutionCoeff)
  info("Using default TV β = 1.0")
  return RegulTV(x0, 1.0)
end

"""
    RegulPM(x₀, T=1.0)

Construct a Perona-Malik regularisation functional with parameter T, and
baseline parameters x₀ (defined in a SolutionCoeff basis).
"""
type RegulPM <: Regul

  ptr::Cxx.CppPtr
  rast::Raster

  function RegulPM(x0::SolutionCoeff, T::Float64)

    rast = x0.rast
    rptr = rast.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new PM(1.0, $(T), &x0, $(rptr));
      return reg;
    """

    regul = new(regulptr, rast)
    finalizer(regul, _regul_delete)
    return regul

  end

end

function RegulPM(x0::SolutionCoeff)
  info("Using default PM T = 1.0")
  return RegulPM(x0, 1.0)
end

# Delete a raster map
_regul_delete{T<:Regul}(regul::T) = finalize(regul.ptr)

# Regul val
"""
  val(regul, x)

Return the value of the regularisation functional defined by `regul` evaluated
with parameters `x` (defined in a SolutionCoeff basis).
"""
function val{T<:Regul}(regul::T, x::SolutionCoeff)

  assert(regul.rast == x.rast)

  xptr = pointer(x.data)
  xlen = length(x)

  icxx"""
    RVector x($xlen, $(xptr), SHALLOW_COPY);
    return $(regul.ptr)->GetValue(x);
  """

end

# grad of regul
"""
  grad(regul, x)

Return the gradient of the regularisation functional defined by `regul` evaluated
with parameters `x` (defined in a SolutionCoeff basis).
"""
function grad{T<:Regul}(regul::T, x::SolutionCoeff)

  assert(regul.rast == x.rast)

  g = SolutionCoeff(regul.rast)

  gptr = pointer(g.data)
  xptr = pointer(x.data)
  xlen = length(x)

  icxx"""
    RVector x($xlen, $xptr, SHALLOW_COPY);
    RVector g($xlen, $gptr, SHALLOW_COPY);
    g = $(regul.ptr)->GetGradient(x);
  """

  return g

end

# Hessian of regul
"""
  val(regul, x)

Return the Hessian of the regularisation functional defined by `regul` evaluated
with parameters `x` (defined in a SolutionCoeff basis).
"""
function hess{T<:Regul}(regul::T, x::SolutionCoeff)
  #
  # assert(regul.rast == x.rast)
  #
  # xptr = pointer(x)
  # xlen = length(x)

  # icxx"""
  #
  #   RVector x(xlen, xvec, SHALLOW_COPY);
  #   RCompRowMatrix *H = new RCompRowMatrix;
  #
  #   int nprm = regul->GetNParam();
  #   int n = x.Dim();
  #   int n0 = n/nprm;
  #   int i, j;
  #
  #   RCompRowMatrix Hi, Hij;
  #   for (i = 0; i < nprm; i++) {
  #       for (j = 0; j < nprm; j++) {
  #           if (!j) {
  #               Hi.New(n0,n0);
  #               if (j==i) regul->SetHess1 (Hi, x, j);
  #           } else {
  #               Hij.New(n0,n0);
  #               if (j==i) regul->SetHess1 (Hij, x, j);
  #               Hi = cath (Hi, Hij);
  #           }
  #       }
  #       if (!i) *H = Hi;
  #       else    *H = catv (*H, Hi);
  #   }
  #
  #   int nval = H->nVal();
  #
  #   lengthptrs[0] = nval;
  #   lengthptrs[1] = n0+1;
  #
  #   return H;
  # """
  #
  # icxx"""
  #   memcpy(colidx, H->colidx, H->nVal()*sizeof(int));
  #   memcpy(rowptr, H->rowptr, (H->nRows()+1)*sizeof(int));
  #   memcpy(val, H->ValPtr(), H->nVal()*sizeof(double));
  #
  #   delete H;
  # """

end
