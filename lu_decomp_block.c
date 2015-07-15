#include "common.h"

static void pfact_inner(int j, const int width, const int N, const int NB, const int NBMIN,
			const int block_width, const int LD, double (*A)[LD], double *b, int *ip)
{
  int local_j  = j % NB;
  //  printf("PANEL : y = %2d - %2d  x = %2d - %2d\n", j, N, j, j + width);

  for(int k=0;k<width;k++,j++,local_j++){
    /* Search pivot * */
    timer_start(PANEL_PIVOT);
    ip[local_j] = (int)cblas_idamax(N-j, &A[j][j], 1) + j;
    timer_stop(PANEL_PIVOT);

    /* SWAP in only Block */
    timer_start(PANEL_SWAP);
    if(ip[local_j] != j)
      cblas_dswap(block_width, &A[j-local_j][j], LD, &A[j-local_j][ip[local_j]], LD);
    timer_stop(PANEL_SWAP);

    /* PANEL Factorization */
    timer_start(PANEL_DSCAL);
    cblas_dscal(N-j-1, 1.0/A[j][j], &A[j][j+1], 1);
    timer_stop(PANEL_DSCAL);

    timer_start(PANEL_DGER);
    cblas_dger(CblasColMajor, N-1-j, width-1-k, -1.0, &A[j][j+1], 1,
	       &A[j+1][j], LD, &A[j+1][j+1], LD);
    timer_stop(PANEL_DGER);
  }
}

static int offset = 0;
static void pdfact(const int j, const int width, const int N, const int NB, const int NBMIN,
		   const int block_width, const int LD, double (*A)[LD], double *b, int *ip)
{
  if(width <= NBMIN){
    pfact_inner(j, width, N, NB, NBMIN, block_width, LD, A, b, ip);
    offset += width;
    return;
  }

  int new_width = width/2;
  int rest      = width - new_width;

  pdfact(j, new_width, N, NB, NBMIN, block_width, LD, A, b, ip);

  timer_start(PANEL_DTRSM);
  //  printf("DTRSM : y = %2d - %2d  x = %2d - %2d\n", j, offset, offset, offset+rest);
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
	      CblasUnit, new_width, rest, 1.0, &A[j][j], LD, &A[offset][j], LD);
  timer_stop(PANEL_DTRSM);

  timer_start(PANEL_DGEMM);
  //  printf("DGEMM : y = %2d - %2d  x = %2d - %2d\n\n", offset, N, offset, offset+rest);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N-offset, rest, new_width,
  	      -1.0, &A[j][offset], LD, &A[offset][j], LD, 1.0, &A[offset][offset], LD);
  timer_stop(PANEL_DGEMM);

  pdfact(offset, rest, N, NB, NBMIN, block_width, LD, A, b, ip);
}

void lu_decomp(const int N, const int NB, const int NBMIN, const int LD, double (*A)[LD], double *b)
{
  timer_start(LU);
  int ip[NB];
  
  for(int j=0;j<N;j+=NB){
    int width = Mmin(j+NB, N) - j;

    timer_start(PANEL);
    pdfact(j, width, N, NB, NBMIN, width, LD, A, b, ip);
    timer_stop(PANEL);

    timer_start(SWAP);
    // SWAP
    for(int k=0;k<width;k++)
      if(j+k != ip[k])
	cblas_dswap(N-j-width+1, &A[j+width][j+k], LD, &A[j+width][ip[k]], LD);
    timer_stop(SWAP);

    // PANEL UPDATE
    //    printf("DTRSM2: y = %2d - %2d  x = %2d - %2d\n", j, j+width, j+width, N+1);
    timer_start(DTRSM);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasUnit, width, N-j-width+1, 1.0, &A[j][j], LD, &A[j+width][j], LD);
    timer_stop(DTRSM);

    timer_start(DGEMM);
    // SUB-MATRIX UPDATE
    //    printf("DGEMM2: y = %2d - %2d  x = %2d - %2d\n\n", j+width, N, j+width, N+1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N-j-width, N-j-width+1, width,
    		-1.0, &A[j][j+width], LD, &A[j+width][j], LD, 1.0, &A[j+width][j+width], LD);
    timer_stop(DGEMM);
  }
  timer_stop(LU);
}
