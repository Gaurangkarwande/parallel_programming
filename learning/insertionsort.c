#include <stdio.h>
#include <omp.h>
#include "utils.h"

#define NUM_THREADS 8

//Cannot be parallelized

void insertionsort_ser(int dim, int arr[])
{
    int i,j, key;
    for (i=1; i < dim; i++)
    {
        key = arr[i];
        j = i-1;
        while (j >= 0 && arr[j] > key)
        {
            arr[j+1] = arr[j];
            j--;
        }
        arr[j+1] = key;
    }
}

int main()
{
    int dim = 8;
    double start, time_s, time_p;
    int* A = malloc(dim * sizeof(int));
    int* B = malloc(dim * sizeof(int));
    initialize_1Darray(dim, A);
    initialize_1Darray(dim, B);

    omp_set_num_threads(NUM_THREADS);
    start = omp_get_wtime();
    insertionsort_ser(dim, B);
    time_s = omp_get_wtime() - start;

    print_1Darray(dim, B);
    
    return 0;
}