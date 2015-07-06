ifeq ($(arch),K)
CC=fccpx
CFLAGS=-Kfast,openmp -SSL2BLAMP -I./COMMON/
CLINK=~/archives/CBLAS/lib/cblas_LINUX.a
else
CC=gcc
CFLAGS=-O2 -std=c99 -I./COMMON/
CLINK=/usr/lib/libcblas.dylib
endif

simple: simple.o common.o main.o hpl_timer.o
	$(CC) $(CFLAGS) $^ $(CLINK) -o $@

blocking: blocking.o common.o main.o hpl_timer.o
	$(CC) $(CFLAGS) $(CLINK) $^ -o $@

simple.o: SERIAL/simple.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

blocking.o: SERIAL/blocking.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

common.o: COMMON/common.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

main.o: SERIAL/main.c COMMON/common.h
	$(CC) $(CFLAGS) $< -c

hpl_timer.o: COMMON/hpl_timer.c COMMON/hpl_timer.h
	$(CC) $(CFLAGS) $< -c

clean:
	rm -f *.o
