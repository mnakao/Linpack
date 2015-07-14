#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/time.h>
#include "hpl_timer.h"
#include "cblas.h"

#define ALIGN 8

#define Mmax(a_,b_) (((a_)>(b_))?(a_):(b_))
#define Mmin(a_,b_) (((a_)<(b_))?(a_):(b_))
#define VERIFY(d_)  (16.0 > d_) ? "PASS"  : "ERROR"
#define MAT() {printf("\n"); print_mat(N, LD, A);}
#define SOL() print_sol(N, b);

extern void verification(const int seed, const int n, const int ld, double (*A)[ld], double *b);
extern void set_value(const int seed, const int n, const int ld, double (*A)[ld], double *b);
extern void print_mat(const int n, const int ld, const double (*A)[ld]);
extern void print_sol(const int, const double*);
extern void set_params(int, char**, int*, int*, int*, int*);
extern void lu_decomp(const int n, const int nb, const int nbmin, const int ld, double A[n+1][ld], double *b);
extern void lu_solve(const int n, const int nb, const int ld, double A[n+1][ld], double *b);
extern void print_performance(int n);
