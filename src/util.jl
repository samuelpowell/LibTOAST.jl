# libtoast.jl: interface to the TOAST++ library
# Copyright (C) 2017 Samuel Powell

"""
    CSRsparse(m,n,rowptr,colidx,values)

Construct a Julia compressed sparse column (CSC) matrix from the row pointer,
column index and values of an ``m \times n`` compressed sparse row (CSR) matrix.
"""
function _CSC(m,n,rowptr,colidx,values)

  valptr = pointer(values)
  _CSC(m,n,rowptr,colidx,valptr)

end

"""
The vector of values can be provided as a pointer.
"""
function _CSC(m,n,rowptr,colidx,values::Ptr)

  # From toastmatlab.cc
  # void CopyMatrix (mxArray **array, const RCompRowMatrix &mat)

  # Array of pointers to row index, column pointer, real pointer
  cscval = Array{Float64}(length(colidx))
  rowidx = Array{Int32}(length(colidx))
  colptr = Array{Int32}(n+1)

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

  # Build sparse with 1-based indexing
  return SparseMatrixCSC(n,m,colptr+1,rowidx+1,cscval)

end
