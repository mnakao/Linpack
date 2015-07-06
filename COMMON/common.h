#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/time.h>
#include "hpl_timer.h"
#include "cblas.h"

#define ALIGN (8)

#define Mmax(a_,b_) (((a_)>(b_))?(a_):(b_))
#define VERIFY(d_)  (16.0 > d_) ? "PASS"  : "ERROR"
#define PRINT_MAT() print_mat(n, ld, A);
#define PRINT_SOL() print_sol(n, b);

extern void verification(int n, int ld, double A[n+1][ld], double *b);
extern void set_value(int n, int ld, double A[n+1][ld], double *b);
extern double calc_norm_2(int n, int ld, double A[n+1][ld]);
extern double calc_norm_1(int n, double *b);
extern void print_mat(int n, int ld, double A[n+1][ld]);
extern void print_sol(int, double*);
extern void set_params(int, char**, int*, int*, int*);
extern void lu_decomp(int n, int nb, int ld, double A[n+1][ld], double *b);
extern void lu_solve(int n, int nb, int ld, double A[n+1][ld], double *b);
extern void print_performance(int n);
