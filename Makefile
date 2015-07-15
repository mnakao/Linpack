ifeq ($(arch),K)
CC=fccpx
CFLAGS=-Kfast,openmp -SSL2BLAMP -I./COMMON/
CLINK=~/archives/CBLAS/lib/cblas_LINUX.a
else
CC=gcc
CFLAGS=-O2 -std=c99 -I./COMMON/
CLINK=/usr/lib/libcblas.dylib
endif

all: block
simple: lu_decomp_simple.o common.o main.o hpl_timer.o lu_solve.o
	$(CC) $(CFLAGS) $^ $(CLINK) -o $@

block: lu_decomp_block.o common.o main.o hpl_timer_blocking.o lu_solve.o
	$(CC) $(CFLAGS) $(CLINK) $^ -o $@

lu_solve.o: SERIAL/lu_solve.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

lu_decomp_simple.o: SERIAL/lu_decomp_simple.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

lu_decomp_block.o: SERIAL/lu_decomp_block.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

common.o: COMMON/common.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

main.o: SERIAL/main.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

hpl_timer.o: COMMON/hpl_timer.c COMMON/hpl_timer.h
	$(CC) $(CFLAGS) $< -c

hpl_timer_blocking.o: COMMON/hpl_timer.c COMMON/hpl_timer.h
	$(CC) $(CFLAGS) $< -c -D_BLOCKING -o hpl_timer_blocking.o

clean:
	rm -f *.o
