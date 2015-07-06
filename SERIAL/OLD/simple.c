#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 16
double A[N][N], b[N];

void print_mat()
{
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      printf("%.2f ", A[i][j]);
    }
    printf("\n");
  }
}

void print_sol()
{
  for(int i=0;i<N;i++)
    printf("%.2f ", b[i]);

  printf("\n");
}

void lu_decomp()
{
  int ip;
  double alpha, tmp, amax;

  for(int k=0;k<N;k++){
    /* Search pivot */
    amax = fabs(A[k][k]);
    ip = k;
    for(int i=k+1;i<N;i++){
      if(fabs(A[i][k]) > amax){
        amax = fabs(A[i][k]);
        ip = i;
      }
    }

    /* SWAP */
    if(ip != k){
      for(int j=k;j<N;j++){
        tmp = A[k][j];
        A[k][j] = A[ip][j];
        A[ip][j] = tmp;
      }

      tmp = b[k];
      b[k] = b[ip];
      b[ip] = tmp;
    }

    /* 前進消去 */
    for(int i=k+1;i<N;i++){
      alpha = -A[i][k]/A[k][k];
      A[i][k] = alpha;
      for(int j=k+1;j<N;j++)
        A[i][j] += alpha * A[k][j];
      
      b[i] += alpha * b[k];
    }
  }
}

void lu_solve()
{
  double tmp;

  /* 後退代入 */
  b[N-1] /= A[N-1][N-1];
  for(int k=N-2;k>=0;k--){
    tmp = b[k];
    for(int j=k+1;j<N;j++)
      tmp -= A[k][j] * b[j];
    b[k] = tmp/A[k][k];
  }
}

int main(int argc, char *argv[])
{
  double sum;

  for(int i=0;i<N;i++){
    sum = 0;
    for(int j=0;j<N;j++){
      A[i][j] = (double)rand()/((double)RAND_MAX+1);
      sum += A[i][j];
    }
    b[i] = sum;
  }

  lu_decomp();
  lu_solve();

  //  print_mat();
  print_sol();

  return 0;
}
