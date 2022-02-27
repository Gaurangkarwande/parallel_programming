#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "utils.h"

#define NUM_THREADS 2

int level = 0;

int partition(int low, int high, int arr[])
{
    int i, j, pivot;
    pivot = arr[high]; i = low - 1;
    for (j = low; j <= high; j++)
    {
        if (arr[j] < pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i+1], &arr[high]);
    return i + 1;
}

int partition_r(int low, int high, int arr[])
{
    srand(213);
    int random = low + rand() % (high - low);
    swap(&arr[random], &arr[high]);
    return partition(low, high, arr);
}

void quicksort_ser(int low, int high, int arr[])
{
    if (low < high)
    {
        int pivot = partition_r(low, high, arr);
        quicksort_ser(low, pivot-1, arr);
        quicksort_ser(pivot+1, high, arr);
    }
}

void quicksort_par(int low, int high, int arr[])
{
    if (low < high)
    {
        int pivot = partition(low, high, arr);
        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task
                {
                    quicksort_par(low, pivot-1, arr);//Thread 1
                }
                #pragma omp task
                {
                    quicksort_par(pivot+1, high, arr);//Thread2
                }
            }
        }
    }
}

// void quicksort(int low, int high, int arr[])
// {
//     printf("%d\n", level++);
//     if (low < high)
//     {
//         int pivot = partition_r(low, high, arr);
//         #pragma omp single
//         {
//             #pragma omp task
//                 quicksort(low, pivot-1, arr);
//             #pragma omp task
//                 quicksort(pivot+1, high, arr);
//         }
//     }
// }

// void sort(int dim, int arr[])
// {
//     #pragma omp parallel shared(dim, arr)
//     {
//         quicksort(0, dim-1, arr);
//     }
// }

int main()
{
    int dim = 0.5e3;
    double start, time_s, time_p;
    int* A = malloc(dim * sizeof(int));
    int* B = malloc(dim * sizeof(int));
    initialize_1Darray(dim, A);
    initialize_1Darray(dim, B);

    omp_set_num_threads(NUM_THREADS);
    start = omp_get_wtime();
    quicksort_ser(0, dim-1, B);
    time_s = omp_get_wtime() - start;

    start = omp_get_wtime();
    quicksort_par(0, dim-1, A);
    time_p = omp_get_wtime() - start;

    printf("time for serial: %f, time for %d threads: %f\n", time_s, NUM_THREADS, time_p);
    
    return 0;
}