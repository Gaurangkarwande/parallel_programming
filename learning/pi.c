#include <stdio.h>
#include <sys/time.h>
#include <omp.h>
#define PAD 16
#define NUM_THREADS 8

static long num_steps = 1000000;
double step;

void calculate_pi(double *pi)
{
    int i, num_threads;
    double sum[NUM_THREADS][PAD];
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        double x;
        int nthreads, thread_ID;
        nthreads = omp_get_num_threads();
        thread_ID = omp_get_thread_num();
        if(thread_ID == 0) num_threads = nthreads;
        for (i = thread_ID, sum[thread_ID][0] = 0.0; i < num_steps; i+= nthreads)
        {
            x = (i+0.5)*step;
            sum[thread_ID][0] += 4.0/(1.0 + x*x);
        }
    }
    for (i = 0; i < num_threads; i++)
        *pi += step*sum[i][0];
}

void calculate_pi_critical(double *pi)
{
    int i;
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        double x;
        int nthreads, thread_ID;
        double thread_sum = 0.0;
        nthreads = omp_get_num_threads();
        thread_ID = omp_get_thread_num();
        for (i = thread_ID; i < num_steps; i+= nthreads)
        {
            x = (i+0.5)*step;
            thread_sum += 4.0/(1.0 + x*x);          //why does directly puttin pi here not work
        }
        printf("Thread id: %d, sum = %f \n", thread_ID, thread_sum);
        #pragma omp critical
            *pi += step*thread_sum;
    }
}

int main()
{
    double time, start, end, pi = 0.0;
    start = omp_get_wtime();
    calculate_pi_critical(&pi);
    end = omp_get_wtime();
    time = end - start;
    printf("Num Threads: %d: pi is %f and Time elapsed: %f ms \n",NUM_THREADS, pi, time);
    return 0;
} 