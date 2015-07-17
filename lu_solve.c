#include "common.h"

//#define _BLAS
#define _BLAS_BLOCK
//#define _OUTER_PRODUCT
//#define _INNER_PRODUCT
void lu_solve(const int N, const int NB, const int LD, double (*A)[LD], double *b)
{
  timer_start(DSOLVE);

#if defined _BLAS
  cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
	      N, &A[0][0], LD, b, 1);
#elif defined _BLAS_BLOCK
  double x[NB];

  int rest = N % NB;
  rest = (rest == 0)? NB : rest;
  int k = N - rest;

  cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
              rest, &A[k][k], LD, &b[k], 1);
  cblas_dcopy(rest, &b[k], 1, x, 1);

  cblas_dgemv(CblasColMajor, CblasNoTrans, NB, rest, -1.0,
	      &A[k][k-NB], LD, x, 1, 1.0, &b[k-NB], 1);
  cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
	      NB, &A[k-NB][k-NB], LD, &b[k-NB], 1);
  cblas_dcopy(NB, &b[k-NB], 1, x, 1);

  cblas_dgemv(CblasColMajor, CblasNoTrans, k-NB, rest, -1.0,
	      &A[k][0], LD, &b[k], 1, 1.0, &b[0], 1);

  for(k=N-rest-NB;k>=NB;k-=NB){
    cblas_dgemv(CblasColMajor, CblasNoTrans, NB, NB, -1.0,
                &A[k][k-NB], LD, x, 1, 1.0, &b[k-NB], 1);
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		NB, &A[k-NB][k-NB], LD, &b[k-NB], 1);
    cblas_dcopy(NB, &b[k-NB], 1, x, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, k-NB, NB, -1.0,
		&A[k][0], LD, &b[k], 1, 1.0, &b[0], 1);
  }

#elif defined _OUTER_PRODUCT
  for(int i=N-1;i>=0;i--){
    double tmp = b[i] /= A[i][i];
    for(int k=0;k<i;k++)
      b[k] -= tmp * A[i][k];

  }
#elif defined _INNER_PRODUCT
  b[N-1] /= A[N-1][N-1];
  for(int k=N-2;k>=0;k--){
    double tmp = b[k];
    for(int j=k+1;j<N;j++)
      tmp -= A[j][k] * b[j];

    b[k] = tmp/A[k][k];
  }
#endif

  timer_stop(DSOLVE);
}
