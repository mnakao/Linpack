#define NUM_TIMERS 13
#define TOTAL       0
#define LU          1
#define DSOLVE      2
#define PANEL       3
#define PANEL_PIVOT 4
#define PANEL_SWAP  5
#define PANEL_DTRSM 6
#define PANEL_DGEMM 7
#define PANEL_DSCAL 8
#define PANEL_DGER  9
#define SWAP        10
#define DTRSM       11
#define DGEMM       12

double start[NUM_TIMERS], elapsed[NUM_TIMERS];

void timer_clear();
void timer_start( int n );
void timer_stop( int n );
void timer_print();
double timer_read( int n );

