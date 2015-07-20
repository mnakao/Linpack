#include "common.h"

void set_params(int argc, char *argv[], int *n, int *nb, int *nbmin, int *ld)
{
  int tmp_n = 16, tmp_ld, tmp_nb = 4, tmp_nbmin = 2;

  switch (argc) {
  case 1:
    break;
  case 2:
    tmp_n = atoi(argv[1]);
    break;
  case 3:
    tmp_n = atoi(argv[1]); tmp_nb = atoi(argv[2]);
    break;
  default:
    tmp_n = atoi(argv[1]); tmp_nb = atoi(argv[2]); tmp_nbmin = atoi(argv[3]);
    break;
  }

  // ld, leading dimension, is a multiple of ALIGN and not a power of 2.
  tmp_ld = ((Mmax(1, tmp_n)-1)/ALIGN) * ALIGN;
  int ip2, i;
  do{
    i = (tmp_ld += ALIGN); ip2 = 1;
    while( i > 1 ) { i >>= 1; ip2 <<= 1; }
  }
  while(tmp_ld == ip2);

  if(tmp_n <= 0 || tmp_nb <= 0 || tmp_nb > tmp_n || tmp_nb < tmp_nbmin){
    printf("Invalid parameter n = %d nb = %d nbmin = %d\n", tmp_n, tmp_nb, tmp_nbmin);
    exit(1);
  }
  else
    printf("n = %d nb = %d nbmin = %d ld = %d\n", tmp_n, tmp_nb, tmp_nbmin, tmp_ld);

  fflush(stdout);
  *n = tmp_n; *ld = tmp_ld; *nb = tmp_nb; *nbmin = tmp_nbmin;
}

// calculate ||A[][]||_oo
static double calc_norm_2(const int N, const int LD, double (*A)[LD])
{
  double tmp_norma[N];
  for(int j=0;j<N;j++){
    tmp_norma[j] = 0.0;
    for(int i=0;i<N;i++){
      tmp_norma[j] += fabs(A[i][j]);
    }
  }

  return tmp_norma[cblas_idamax(N, tmp_norma, 1)];
}

// calculate ||b[]||_oo
static double calc_norm_1(const int N, double *b)
{
  return fabs(b[cblas_idamax(N, b, 1)]);
}

// ||Ax-b||_oo/(eps*(||A||_oo*||x||_oo+||b||_oo)*N) < 16.0 ?
void verification(const int seed, const int N, const int LD, double (*A)[LD], double *b)
{
  double x[N], norma, normb, normx;

  memcpy(x, b, sizeof(double)*N);  // copy b to x
  normx = calc_norm_1(N, x);       // ||x||_oo
  set_value(seed, N, LD, A, b);    // reset A and b
  norma = calc_norm_2(N, LD, A);   // ||A||_oo
  normb = calc_norm_1(N, b);       // ||b||_oo

  // b = A*x - b
  cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, &A[0][0], LD, x, 1, -1.0, b, 1);

  double resid = calc_norm_1(N, b);  // ||Ax-b||_oo
  double v = resid/(DBL_EPSILON*(norma * normx + normb)*N);

  printf("||Ax-b||_oo/(eps*(||A||_oo*||x||_oo+||b||_oo)*N) = %.3f (%s)\n", v, VERIFY(v));
  printf("---\n");
  fflush(stdout);
}

// Initialize A[][] & b[]
void set_value(const int seed, const int N, const int LD, double (*A)[LD], double *b)
{
  srand48((long int)seed);

  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      A[i][j] = -1 * drand48() + 0.5;

  for(int i=0;i<N;i++)
    b[i] = -1 * drand48() + 0.5;
}

void print_mat(const int N, const int LD, const double (*A)[LD])
{
  for(int j=0;j<N;j++){
    for(int i=0;i<N;i++){
      printf("%.3f ", A[i][j]);
    }
    printf("\n");
  }
}

void print_sol(const int N, const double *b)
{
  for(int i=0;i<N;i++)
    printf("%.3f ", b[i]);

  printf("\n");
}

void print_performance(const int N)
{
  double Gflops = ( ( (double)(N) / 1.0e+9 ) * ( (double)(N) / timer_read(TOTAL) ) ) *
    ( (2.0/3.0) * (double)(N) + (3.0/2.0) );

  printf("%.2f GFlops\n", Gflops);
  printf("---\n");
  fflush(stdout);
}


/* The function HPL_dlaswp00N() is a part of hpl-2.1 in the following URL. */
/* http://www.netlib.org/benchmark/hpl/HPL_dlaswp00N.html */
#define HPL_LASWP00N_DEPTH      32
#define HPL_LASWP00N_LOG2_DEPTH 5
void HPL_dlaswp00N(const int M, const int N, double *A, const int LDA, const int *IPIV )
{
  register double r;
  double          *a0, *a1;
  const int       incA = (int)( (unsigned int)(LDA) <<
					   HPL_LASWP00N_LOG2_DEPTH );
  int             ip, nr, nu;
  register int    i, j;

  if( ( M <= 0 ) || ( N <= 0 ) ) return;

  nr = N - ( nu = (int)( ( (unsigned int)(N) >> HPL_LASWP00N_LOG2_DEPTH )
			 << HPL_LASWP00N_LOG2_DEPTH ) );

  for( j = 0; j < nu; j += HPL_LASWP00N_DEPTH, A += incA ){
    for( i = 0; i < M; i++ ){
      if( i != ( ip = IPIV[i] ) ){
	a0 = A + i; a1 = A + ip;

	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#if ( HPL_LASWP00N_DEPTH >  1 )
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#endif
#if ( HPL_LASWP00N_DEPTH >  2 )
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#endif
#if ( HPL_LASWP00N_DEPTH >  4 )
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#endif
#if ( HPL_LASWP00N_DEPTH >  8 )
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#endif
#if ( HPL_LASWP00N_DEPTH > 16 )
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
	r = *a0; *a0 = *a1; *a1 = r; a0 += LDA; a1 += LDA;
#endif
      }
    }
  }

  if( nr > 0 ){
    for( i = 0; i < M; i++ ){
      if( i != ( ip = IPIV[i]) ){
	a0 = A + i; a1 = A + ip;
	for( j = 0; j < nr; j++, a0 += LDA, a1 += LDA )
	  { r = *a0; *a0 = *a1; *a1 = r; }
      }
    }
  }
}
