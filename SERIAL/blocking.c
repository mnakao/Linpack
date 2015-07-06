#include "common.h"

void lu_decomp(int n, int nb, int ld, double A[n+1][ld], double *b)
{
  for(int k=0;k<n;k++){
    int ip = (int)cblas_idamax(n-k, &A[k][k], 1) + k;       /* Search pivot */

    if(ip != k)
      cblas_dswap(n-k+1, &A[k][k], ld, &A[k][ip], ld);      /* SWAP */
    
    /* UPDATE */
    cblas_dscal(n-k-1, -1/A[k][k], &A[k][k+1], 1);
    cblas_dger(CblasRowMajor, n-k, n-k-1, 1.0, &A[k+1][k], ld,
    	       &A[k][k+1], 1, &A[k+1][k+1], ld);
    //    cblas_dger(CblasColMajor, n-k-1, n-k, 1.0, &A[k][k+1], 1,
    //	       &A[k+1][k], ld, &A[k+1][k+1], ld);
  }
}

void lu_solve(int n, int nb, int ld, double A[n+1][ld], double *b)
{
  double tmp;

  b[n-1] /= A[n-1][n-1];
  for(int k=n-2;k>=0;k--){
    tmp = b[k];
    for(int j=k+1;j<n;j++)
      tmp -= A[j][k] * b[j];

    b[k] = tmp/A[k][k];
  }
}
