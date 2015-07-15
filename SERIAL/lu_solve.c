#include "common.h"

//#define _OUTER_PRODUCT
//#define _OUTER__PRODUCT_BLOCKING
//#define _INNER_PRODUCT
//#define _INNER_PRODUCT_BLOCKING
void lu_solve(const int N, const int NB, const int LD, double (*A)[LD], double *b)
{
  timer_start(DSOLVE);

  cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
	      N, &A[0][0], LD, b, 1);

#if defined _OUTER_PRODUCT
  for(int i=N-1;i>=0;i--){
    double tmp = b[i] /= A[i][i];
    for(int k=0;k<i;k++)
      b[k] -= tmp * A[i][k];
  }
#elif defined _OUTER__PRODUCT_BLOCKING
  for(int i=N-1;i>=0;i-=NB){
    int end = Mmax(i-NB+1, 0);
    for(int j=i;j>=end;j--){
      double tmp = b[j] /= A[j][j];
      for(int k=0;k<j;k++)
	b[k] -= tmp * A[j][k];
    }
  }
#elif defined _INNER_PRODUCT
  b[N-1] /= A[N-1][N-1];
  for(int k=N-2;k>=0;k--){
    double tmp = b[k];
    for(int j=k+1;j<N;j++){
      tmp -= A[j][k] * b[j];
    }
    b[k] = tmp/A[k][k];
  }
#elif defined _INNER_PRODUCT_BLOCKING
  b[N-1] /= A[N-1][N-1];
  for(int k=N-2;k>=0;k-=NB){
    int end = Mmax(k-NB+1, 0);
    for(int i=k;i>=end;i--){
      double tmp = b[i];
      for(int j=i+1;j<N;j++){
        tmp -= A[j][i] * b[j];
      }
      b[i] = tmp/A[i][i];
    }
  }
#endif

  timer_stop(DSOLVE);
}
