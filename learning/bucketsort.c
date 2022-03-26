#include <stdio.h>
#include <omp.h>
#include "utils.h"
#include "quicksort.h"

#define NUM_THREADS 8

void initialize_matrix(int row, int col, int mat[])
{
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            mat[i*col + j] = -1;
        }
    }
   return; 
}

int compare_array(int dim, int first[], int second[])
{
    for (int i = 0; i < dim; ++i)
        if (first[i] != second[i])
            return -1;
    return 0;
}

void print_buckets(int n_buckets, int dim, int buckets[], int bucket_sizes[])
{
    int i, j;
    printf("\n Buckets are: \n");
    for ( i = 0; i < n_buckets; i++) {
        for (j = 0; j < bucket_sizes[i]; j++) {
            printf("%d \t", buckets[i*dim + j]);
        }
        printf("\n");
    }
   return; 
}

void initialize_array(int dim, int mat[])
{
    int i;
    for ( i = 0; i < dim; i++) {
        mat[i] = 0;
    }
   return; 
}

void insert_bucket(int num, int bucketID, int dim, int buckets[])
{
    int i = 0;
    while(buckets[bucketID*dim + i] != -1)
        i++;
    buckets[bucketID*dim + i] = num;
}

int get_bucket_size(int bucketID, int dim, int buckets[])
{
    int i = 0;
    while(buckets[bucketID*dim + i] != -1)
        i++;
    return i;
}

void extract_bucket(int id, int bucketID, int dim, int buckets[], int bucket_sizes[], int arr[])
{
    for (int i=0; i < bucket_sizes[bucketID]; i++)
    {
        arr[id] = buckets[bucketID*dim + i];
        id++;
    }
}

void bucketsort_ser(int n_buckets, int dim, int arr[])
{
    int i, bucketID;
    int* buckets = malloc((n_buckets * dim) * sizeof(int));
    int* bucket_sizes = malloc(n_buckets * sizeof(int));
    initialize_matrix(n_buckets, dim, buckets);
    initialize_array(n_buckets, bucket_sizes);
    for (i = 0; i < dim; i++)
    {
        bucketID = arr[i]/ (dim/n_buckets);
        insert_bucket(arr[i], bucketID, dim, buckets);
    }

    for (bucketID=0; bucketID < n_buckets; bucketID++)
    {
        bucket_sizes[bucketID] = get_bucket_size(bucketID, dim, buckets);
    }

    for (bucketID=0; bucketID < n_buckets; bucketID++)
    {
        quicksort_ser(bucketID*dim, bucketID*dim + bucket_sizes[bucketID] - 1, buckets);
    }

    i = 0;
    for (bucketID=0; bucketID < n_buckets; bucketID++)
    {
        extract_bucket(i, bucketID, dim, buckets, bucket_sizes, arr);
        i += bucket_sizes[bucketID];
    }
    free(buckets); free(bucket_sizes);
}

void bucketsort_par(int n_buckets, int dim, int arr[])
{
    int i, bucketID;
    int* buckets = malloc((n_buckets * dim) * sizeof(int));
    int* bucket_sizes = malloc(n_buckets * sizeof(int));
    initialize_matrix(n_buckets, dim, buckets);
    initialize_array(n_buckets, bucket_sizes);

    #pragma omp parallel shared(buckets, bucket_sizes)
    {
        #pragma omp for private(i, bucketID)
            for (i = 0; i < dim; i++)
            {
                bucketID = arr[i]/ (dim/n_buckets);
                #pragma omp critical
                insert_bucket(arr[i], bucketID, dim, buckets);
            }
        
        #pragma omp for private(bucketID)
            for (bucketID=0; bucketID < n_buckets; bucketID++)
            {
                bucket_sizes[bucketID] = get_bucket_size(bucketID, dim, buckets);
            }

        #pragma omp for private(bucketID)
            for (bucketID=0; bucketID < n_buckets; bucketID++)
            {
                quicksort_ser(bucketID*dim, bucketID*dim + bucket_sizes[bucketID] - 1, buckets);
            }
    }
        
    i = 0;
    for (bucketID=0; bucketID < n_buckets; bucketID++)
    {
        extract_bucket(i, bucketID, dim, buckets, bucket_sizes, arr);
        i += bucket_sizes[bucketID];
    }
    free(buckets); free(bucket_sizes);
}

int main()
{
    int dim, range_low, range_high, n_buckets;
    dim = 1e5 ; range_low = 0, range_high = dim-1;
    n_buckets = NUM_THREADS;
    double start, time_s, time_p;
    int* A = malloc(dim * sizeof(int));
    int* B = malloc(dim * sizeof(int));
    initialize_1Darray_inRange(dim, A, range_low, range_high);
    initialize_1Darray_inRange(dim, B, range_low, range_high);

    omp_set_num_threads(NUM_THREADS);
    start = omp_get_wtime();
    bucketsort_ser(n_buckets, dim, B);
    time_s = omp_get_wtime() - start;

    start = omp_get_wtime();
    bucketsort_par(n_buckets, dim, A);
    time_p = omp_get_wtime() - start;

    if (compare_array(dim, A, B) == -1)
    {
        printf("Wrong sorting. \n");
    }

    printf("time for serial: %f, time for %d threads: %f\n", time_s, NUM_THREADS, time_p);
}