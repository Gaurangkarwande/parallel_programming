#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int r, i;
} complex;

complex complex_multiply(complex a, complex b)
{
    complex c = {a.r*b.r - a.i*b.i, a.r*b.i + a.i*b.r};
    return c;
}

void complex_increment(complex *sum, complex a)
{
    sum->r += a.r;
    sum->i += a.i;
    return;
}

void initialize_matrix(int dim, complex mat[], int result)
{
    int i, j;
    for ( i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            mat[i*dim + j].r = 0;
            mat[i*dim + j].i = 0;
            if (i == j && !result)
            {
                mat[i*dim + j].r = 1;
                mat[i*dim + j].i = 1;
            }
        }
    }
   return; 
}

void transpose_matrix(int dim, complex mat[], complex result[])
{
    int i, j;
    for ( i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            result[j*dim + i] = mat[i*dim + j];
   return; 
}

void print_matrix(int row, int col, complex mat[])
{
    int i, j;
    printf("\n");
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d + j%d \t", mat[i*col + j].r, mat[i*col + j].i);
        }
        printf("\n");
    }
    printf("\n");
   return; 
}

int compare_matrix(int dim, complex first[], complex second[])
{
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) 
            if (first[i*dim + j].r != second[i*dim + j].r || first[i*dim + j].i != second[i*dim + j].i)
                return -1;
    return 0;
}

void multiply_naive(complex first[], complex second[], complex result[], int dim)
{
    int i,j,k;
    complex sum;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            sum.r = 0;
            sum.i = 0;
            for (k = 0; k < dim; k++) {
                complex_increment(&sum, complex_multiply(first[i*dim + k], second[k*dim + j]));
            }
            result[i*dim + j] = sum;
        }
    }
}

void multiply_transpose(complex first[], complex second[], complex result[], int block_dim, int dim)
{
    //Second matrix is transposed

    int i,j,k;
    complex sum;
    for (i = 0; i < block_dim; i++) {
        for (j = 0; j < block_dim; j++) {
            sum.r = 0;
            sum.i = 0;
            for (k = 0; k < dim; k++) {
                complex_increment(&sum, complex_multiply( first[i*dim + k], second[j*dim + k]));
            }
            result[i*block_dim + j] = sum;
        }
    }
}

void multiply_submatrix(complex first[], complex second[], complex result[], int block_dim, int dim)
{
    //Second matrix is transposed

    int i,j,k;
    complex sum;
    for (i = 0; i < block_dim; i++) {
        for (j = 0; j < block_dim; j++) {
            sum.r = 0;
            sum.i = 0;
            for (k = 0; k < dim; k++) {
                complex_increment(&sum, complex_multiply( first[i*dim + k], second[j*dim + k]));
            }
            result[i*dim + j] = sum;
        }
    }
}

int main(int argc, char *argv[])
{
    int dim = 840;
    double start, time_s;
    complex* A = malloc((dim * dim) * sizeof(complex));
    complex* B = malloc((dim * dim) * sizeof(complex));
    complex* Baseline = malloc((dim * dim) * sizeof(complex));

    initialize_matrix(dim, A, 0);
    initialize_matrix(dim, B, 0);
    initialize_matrix(dim, Baseline, 1);
    
    start = MPI_Wtime();
    multiply_naive(A, B, Baseline, dim);
    time_s = MPI_Wtime() - start;
    printf("For matrix size %d : Time for serial: %f ms \n", dim, time_s*1000);

    return 0;
}