#include "common.h"

int main(int argc, char *argv[])
{
  double *vptr = NULL;
  int n, ld, nb, nbmin, seed = 1;

  set_params(argc, argv, &n, &nb, &nbmin, &ld);

  /* malloc A[n+1][ld] */
  vptr = (double *)malloc(sizeof(double) * ((n+1) * ld + ALIGN));
  double (*A)[ld] = (double (*)[ld])HPL_PTR( vptr, ((size_t)(8) * sizeof(double) ) );

  /* A[n][0..n-1] is x[0..n-1] */
  double *b = &A[n][0];

  set_value(seed, n, ld, A, b);

  timer_clear();
  timer_start(TOTAL);
  lu_decomp(n, nb, nbmin, ld, A, b);
  lu_solve(n, nb, ld, A, b);
  timer_stop(TOTAL);

  timer_print();
  print_performance(n);

  verification(seed, n, ld, A, b);

  free(vptr);
  return 0;
}
