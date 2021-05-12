/* mountain.c - Generate the memory mountain. */
/* $begin mountainmain */
#include <stdlib.h>
#include <stdio.h>
#include "fcyc2.h" /* K-best measurement timing routines */
#include "clock.h" /* routines to access the cycle counter */

#define MINBYTES (1 << 1) /* Working set size ranges from 1 KB */
#define MAXBYTES (1 << 15) 
#define MAXSTRIDE 32       /* Strides range from 1 to 32 */
#define STRIDESTRIDE 2  /* increment stride by this amount each time */
#define MAXELEMS MAXBYTES / sizeof(long)

long A[MAXELEMS][MAXELEMS]; /* The array we'll be traversing */
/* $end mountainmain */
void init_data(long A[][MAXELEMS], int n);
int test(int elems, int stride);
double run(int size, int stride, double Mhz);

/* $begin mountainmain */
int main()
{
    int size;   /* Working set size (in bytes) */
    int stride; /* Stride (in array elements) */
    double Mhz; /* Clock frequency */
    init_data(A, MAXELEMS); /* Initialize each element in data to 1 */
    Mhz = mhz(0);              /* Estimate the clock frequency */
    /* $end mountainmain */
    /* Not shown in the text */
    printf("Clock frequency is approx. %.1f MHz\n", Mhz);
    printf("Memory mountain (MB/sec)\n");

    printf("\t");
    for (stride = 1; stride <= MAXSTRIDE; stride += STRIDESTRIDE)
        printf("s%d\t", stride);
    printf("\n");

    /* $begin mountainmain */
    for (size = MAXBYTES; size >= MINBYTES; size >>= 1)
    {
        /* $end mountainmain */
        /* Not shown in the text */
        if (size > (1 << 20))
            printf("%dm\t", size / (1 << 20));
        else
            printf("%dk\t", size / 1024);

        /* $begin mountainmain */
        for (stride = 1; stride <= MAXSTRIDE; stride += STRIDESTRIDE)
        {
            printf("%.0f\t", run(size, stride, Mhz));
        }
        printf("\n");
    }
    exit(0);
}
/* $end mountainmain */

/* init_data - initializes the array */
void init_data(long A[][MAXELEMS],  int n)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = 1;
        }
    }
}

/* $begin mountainfuns */
int test(int elems, int stride) /* The test function */
{
    int i, k, j, r, result = 0;
    volatile int sink;

    for (k = 0; k < elems; k+= stride)
    {
        for (i = 0; i < elems; i+= stride)
        {
            r = A[i][k];
            for (j = 0; j < elems; j+= stride)
                result += r * A[k][j];
        }
    }


    sink = result; /* So compiler doesn't optimize away the loop */
    return sink;
}

/* Run test(elems, stride) and return read throughput (MB/s) */
double run(int size, int stride, double Mhz)
{
    double cycles;
    int elems = size / sizeof(long);

    test(elems, stride);                     /* warm up the cache */
    cycles = fcyc2(test, elems, stride, 0);  /* call test(elems,stride) */
    return (size / stride) / (cycles / Mhz); /* convert cycles to MB/s */
}
/* $end mountainfuns */
