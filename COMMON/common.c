#include "common.h"

void set_params(int argc, char *argv[], int *n, int *nb, int *ld)
{
  int tmp_n, tmp_ld, tmp_nb;

  if(argc == 1){ tmp_n = 16; tmp_nb = 4; }
  else if(argc == 2) { tmp_n = atoi(argv[1]); tmp_nb = 4; }
  else { tmp_n = atoi(argv[1]); tmp_nb = atoi(argv[2]); }

  // ld, leading dimension, is a multiple of ALIGN and not a power of 2.
  tmp_ld = ((Mmax(1, tmp_n)-1)/ALIGN) * ALIGN;
  int ip2, i;
  do{
    i = (tmp_ld += ALIGN); ip2 = 1;
    while( i > 1 ) { i >>= 1; ip2 <<= 1; }
  }
  while(tmp_ld == ip2);

  if(tmp_n <= 0 || tmp_nb <= 0){
    printf("Invalid parameter n = %d nb = %d\n", tmp_n, tmp_nb);
    exit(1);
  }
  else
    printf("n = %d nb = %d ld = %d\n", tmp_n, tmp_nb, tmp_ld);

  fflush(stdout);
  *n = tmp_n; *ld = tmp_ld; *nb = tmp_nb;
}

// ||Ax-b||_oo/(eps*(||A||_oo*||x||_oo+||b||_oo)*N) < 16.0 ?
void verification(int n, int ld, double A[n+1][ld], double *b)
{
  double x[n], tmp, norma, normb, normx;

  memcpy(x, b, sizeof(double)*n);  // copy b to x
  normx = calc_norm_1(n, x);       // ||x||_oo
  set_value(n, ld, A, b);          // reset A and b
  norma = calc_norm_2(n, ld, A);   // ||A||_oo
  normb = calc_norm_1(n, b);       // ||b||_oo

  // b = A*x - b
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, &A[0][0], ld, x, 1, -1.0, b, 1);

  double resid = calc_norm_1(n, b);  // ||Ax-b||_oo
  double v = resid/(DBL_EPSILON*(norma * normx + normb)*n);

  printf("||Ax-b||_oo/(eps*(||A||_oo*||x||_oo+||b||_oo)*N) = %.3f (%s)\n", v, VERIFY(v));
}

// Initialize A[][] & b[]
void set_value(int n, int ld, double A[n+1][ld], double *b)
{
  srand48((long int)1);

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      A[i][j] = -1 * drand48() + 0.5;

  for(int i=0;i<n;i++)
    b[i] = -1 * drand48() + 0.5;
}

// calculate ||A[][]||_oo
double calc_norm_2(int n, int ld, double A[n+1][ld])
{
  double tmp_norma[n];
  for(int j=0;j<n;j++){
    tmp_norma[j] = 0.0;
    for(int i=0;i<n;i++){
      tmp_norma[j] += fabs(A[i][j]);
    }
  }

  return tmp_norma[cblas_idamax(n, tmp_norma, 1)];
}

// calculate ||b[]||_oo
double calc_norm_1(int n, double *b)
{
  return fabs(b[cblas_idamax(n, b, 1)]);
}

void print_mat(int n, int ld, double A[n+1][ld])
{
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%.3f ", A[i][j]);
    }
    printf("\n");
  }
}

void print_sol(int n, double *b)
{
  for(int i=0;i<n;i++)
    printf("%.3f ", b[i]);

  printf("\n");
}

void print_performance(int n)
{
  double Gflops = ( ( (double)(n) / 1.0e+9 ) * ( (double)(n) / timer_read(TOTAL) ) ) *
    ( (2.0/3.0) * (double)(n) + (3.0/2.0) );
  printf("%.2f GFlops\n", Gflops);
}
