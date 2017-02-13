# TOAST.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

# Imports

# Export types
export Regul, RegulTK0, RegulTK1, RegulTV, RegulPM

# Export methods
export val, grad, hess

abstract Regul

type RegulTK0 <: Regul

  ptr::Cxx.CppPtr
  rmap::RasterMap

  function RegulTK0(x0::SolutionCoeff)

    rmap = x0.rmap
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      RVector xs($(xlen), 1);
      Regularisation *reg = new Tikhonov0 (1.0, &x0, &xs);
      return reg;
    """

    regul = new(regulptr, rmap)
    finalizer(regul, _regul_delete)
    return regul
  end

end

type RegulTK1 <: Regul

  ptr::Cxx.CppPtr
  rmap::RasterMap

  function RegulTK1(x0::SolutionCoeff)

    rmap = x0.rmap
    rptr = rmap.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      void *kapref = 0;
      bool istensor = false;
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new TK1 (1.0, &x0, $(rptr), kapref, istensor);
      return reg;
    """

    regul = new(regulptr, rmap)
    finalizer(regul, _regul_delete)
    return regul
  end

end

type RegulTV <: Regul

  ptr::Cxx.CppPtr
  rmap::RasterMap

  function RegulTV(x0::SolutionCoeff, β::Float64)

    rmap = x0.rmap
    rptr = rmap.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new TV(1.0, $(β), &x0, $(rptr));
      return reg;
    """

    regul = new(regulptr, rmap)
    finalizer(regul, _regul_delete)
    return regul

  end

end

function RegulTV(x0::SolutionCoeff)
  info("Using default TV β = 1.0")
  return RegulTV(x0, 1.0)
end


type RegulPM <: Regul

  ptr::Cxx.CppPtr
  rmap::RasterMap

  function RegulPM(x0::SolutionCoeff, T::Float64)

    rmap = x0.rmap
    rptr = rmap.ptr
    xptr = pointer(x0.data)
    xlen = length(x0)

    regulptr = icxx"""
      RVector x0($(xlen), $(xptr), SHALLOW_COPY);
      Regularisation *reg = new PM(1.0, $(T), &x0, $(rptr));
      return reg;
    """

    regul = new(regulptr, rmap)
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
function val{T<:Regul}(regul::T, x::SolutionCoeff)

  assert(regul.rmap == x.rmap)

  xptr = pointer(x.data)
  xlen = length(x)

  icxx"""
    RVector x($xlen, $(xptr), SHALLOW_COPY);
    return $(regul.ptr)->GetValue(x);
  """

end

# grad of regul
function grad{T<:Regul}(regul::T, x::SolutionCoeff)

  assert(regul.rmap == x.rmap)

  g = SolutionCoeff(regul.rmap)

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
function hess{T<:Regul}(regul::T, x::SolutionCoeff)
  #
  # assert(regul.rmap == x.rmap)
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
