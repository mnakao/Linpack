#include "common.h"

void lu_decomp(int n, int nb, int ld, double A[n+1][ld], double *b)  /* nb is not used. */
{
  timer_start(LU);

  for(int k=0;k<n;k++){
    timer_start(PIVOT);
    int ip = (int)cblas_idamax(n-k, &A[k][k], 1) + k;       /* Search pivot */
    timer_stop(PIVOT);

    timer_start(SWAP);
    if(ip != k)
      cblas_dswap(n-k+1, &A[k][k], ld, &A[k][ip], ld);      /* SWAP */
    timer_stop(SWAP);
    
    /* PANEL Factorization */
    timer_start(PANEL);
    cblas_dscal(n-k-1, -1/A[k][k], &A[k][k+1], 1);
    timer_stop(PANEL);

    /* UPDATE */
    timer_start(UPDATE);
    cblas_dger(CblasRowMajor, n-k, n-k-1, 1.0, &A[k+1][k], ld,
    	       &A[k][k+1], 1, &A[k+1][k+1], ld);
    timer_stop(UPDATE);
    //    cblas_dger(CblasColMajor, n-k-1, n-k, 1.0, &A[k][k+1], 1,
    //	       &A[k+1][k], ld, &A[k+1][k+1], ld);
  }
  timer_stop(LU);
}

void lu_solve(int n, int nb, int ld, double A[n+1][ld], double *b)  /* nb is not used. */
{
  timer_start(DSOLVE);

  b[n-1] /= A[n-1][n-1];
  for(int k=n-2;k>=0;k--){
    double tmp = b[k];
    for(int j=k+1;j<n;j++)
      tmp -= A[j][k] * b[j];

    b[k] = tmp/A[k][k];
  }

  timer_stop(DSOLVE);
}
