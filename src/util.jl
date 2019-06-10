# LibTOAST.jl: interface to the TOAST++ library
# Copyright (C) 2019 Samuel Powell

# Construct a one-based Julia CSC matrix from the row pointer, column index
# and values of an m x n zero-based CSR matrix.
function _CSR0_to_CSC1(m,n,rowptr,colidx,values)

  _CSR0_to_CSC1(m,n,rowptr,colidx,pointer(values))
  
end

# Construct a one-based Julia CSC matrix from the row pointer, column index
# and values of an m x n zero-based CSR matrix.
function _CSR0_to_CSC1(m,n,rowptr,colidx,values::Ptr)

  # From toastmatlab.cc
  # void CopyMatrix (mxArray **array, const RCompRowMatrix &mat)

  # Array of pointers to row index, column pointer, real pointer
  cscval = Array{Float64}(undef, length(colidx))
  rowidx = Array{Cint}(undef, length(colidx))
  colptr = Array{Cint}(undef, n+1)

  icxx"""
    int i, j, k;
    int m = $m;
    int n = $n;
    int *rowptr = $(pointer(rowptr));
    int *colidx = $(pointer(colidx));


    const double *pval = $values;
    int nz = rowptr[m];

    // This should maybe the larger than int!
    int *rcount = new int[n];
    int idx;

    //mxArray *tmp = mxCreateSparse (m, n, nz, mxREAL);
    for (i = 0; i < n; i++) rcount[i] = 0;
    for (i = 0; i < nz; i++) rcount[colidx[i]]++;

    double *pr = $(pointer(cscval));
    int *ir = $(pointer(rowidx));
    int *jc = $(pointer(colptr));

    jc[0] = 0;
    for (i = 0; i < n; i++) jc[i+1] = jc[i]+rcount[i];
    for (i = 0; i < n; i++) rcount[i] = 0;
    for (i = 0; i < m; i++)
    for (k = rowptr[i]; k < rowptr[i+1]; k++) {
      j = colidx[k];
      idx = jc[j]+rcount[j];
      ir[idx] = i;
      pr[idx] = pval[k];
      rcount[j]++;
    }

    delete []rcount;
  """

  # Build sparse with 1-based indexing and Int size pointers
  return SparseMatrixCSC(n,m,Vector{Int}(colptr.+1),Vector{Int}(rowidx.+1),cscval)

end
