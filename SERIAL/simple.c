#include "common.h"

void lu_decomp(const int N, const int NB, const int NBMIN, const int LD, double (*A)[LD], double *b)  /* NB is not used. */
{
  timer_start(LU);

  for(int k=0;k<N;k++){
    timer_start(PANEL_PIVOT);
    int ip = (int)cblas_idamax(N-k, &A[k][k], 1) + k;       /* Search pivot */
    timer_stop(PANEL_PIVOT);

    timer_start(SWAP);
    if(ip != k)
      cblas_dswap(N-k+1, &A[k][k], LD, &A[k][ip], LD);      /* SWAP */
    timer_stop(SWAP);

    /* PANEL Factorization */
    timer_start(PANEL);
    cblas_dscal(N-k-1, 1.0/A[k][k], &A[k][k+1], 1);
    timer_stop(PANEL);

    /* UPDATE */
    timer_start(PANEL_DGER);
    cblas_dger(CblasRowMajor, N-k, N-k-1, -1.0, &A[k+1][k], LD,
    	       &A[k][k+1], 1, &A[k+1][k+1], LD);
    timer_stop(PANEL_DGER);
  }
  timer_stop(LU);
}

void lu_solve(const int N, const int NB, const int LD, double (*A)[LD], double *b)  /* NB is not used. */
{
  timer_start(DSOLVE);

  b[N-1] /= A[N-1][N-1];
  for(int k=N-2;k>=0;k--){
    double tmp = b[k];
    for(int j=k+1;j<N;j++){
      tmp -= A[j][k] * b[j];
    }

    b[k] = tmp/A[k][k];
  }

  timer_stop(DSOLVE);
}
