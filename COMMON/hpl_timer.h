#define NUM_TIMERS 7
#define TOTAL  0
#define LU     1
#define DSOLVE 2
#define PIVOT  3
#define SWAP   4
#define PANEL  5
#define UPDATE 6

double start[NUM_TIMERS], elapsed[NUM_TIMERS];

void timer_clear();
void timer_start( int n );
void timer_stop( int n );
void timer_print();
double timer_read( int n );

