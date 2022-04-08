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
            printf("%d + j%d   ", mat[i*col + j].r, mat[i*col + j].i);
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

void multiply_naive(int r1, int c1, complex first[], int r2, int c2, complex second[], complex result[])
{
    int i,j,k;
    complex sum;
    for (i = 0; i < r1; i++) {
        for (j = 0; j < c2; j++) {
            sum.r = 0;
            sum.i = 0;
            for (k = 0; k < c1; k++) {
                complex_increment(&sum, complex_multiply(first[i*c1 + k], second[k*c2 + j]));
            }
            result[i*c2 + j] = sum;
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
    int dim = 1680;
    double start, time_s, time_p;
    complex* A = malloc((dim * dim) * sizeof(complex));
    complex* B = malloc((dim * dim) * sizeof(complex));
    complex* C = malloc((dim * dim) * sizeof(complex));
    complex* Baseline = malloc((dim * dim) * sizeof(complex));
    complex* Bt = malloc((dim * dim) * sizeof(complex));

    initialize_matrix(dim, A, 0);
    initialize_matrix(dim, B, 0);
    initialize_matrix(dim, C, 1 );
    initialize_matrix(dim, Baseline, 1);
    transpose_matrix(dim, B, Bt);
    
    start = MPI_Wtime();
    // multiply_naive(dim, dim, A, dim, dim, B, Baseline);
    time_s = MPI_Wtime() - start;

    start = MPI_Wtime();
    complex *sub_A, *sub_B, *sub_C;
    int rank, numprocs, block_dim, b_rows, b_cols, i, j, k, proc_i, proc_j;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    MPI_Datatype c_type;
    MPI_Type_contiguous(2, MPI_INT, &c_type);
    MPI_Type_commit(&c_type);

    block_dim = dim/sqrt(numprocs);
    b_rows = b_cols = (int) sqrt(numprocs);

    if (rank == 0)  
    {
        MPI_Datatype row_vec;
        MPI_Type_vector(dim, 1, 1, c_type, &row_vec);
        MPI_Type_commit(&row_vec);

        MPI_Datatype col_vec;
        MPI_Type_vector(dim, block_dim, dim, c_type, &col_vec);
        MPI_Type_commit(&col_vec);

        MPI_Request request_A[numprocs-1];
        MPI_Request request_B[numprocs-1];

        sub_C = malloc((block_dim*block_dim) * sizeof(complex));
        // print_matrix(dim, dim, C);
        for (i = 0; i < b_rows; i++)
        {
            for (j = 0; j < b_cols; j++)
            {
                if (i == 0 && j == 0)
                    continue;
                MPI_Isend(A+i*block_dim*dim, block_dim, row_vec, i*b_rows+j, 0, MPI_COMM_WORLD, &request_A[ i*b_rows+j - 1]);
                MPI_Isend(B+j*block_dim, 1, col_vec, i*b_rows+j, 1, MPI_COMM_WORLD, &request_B[ i*b_rows+j - 1]);
                
            }
            
        }
        multiply_submatrix(A, Bt, C, block_dim, dim);
        // print_matrix(dim, dim, C);
        for (k = 1; k < numprocs; k++)
        {
            MPI_Wait(&request_A[k-1], MPI_STATUS_IGNORE);
            MPI_Wait(&request_B[k-1], MPI_STATUS_IGNORE);
            proc_i = (int) k/b_rows;
            proc_j = k%b_cols;
            if (proc_i == 0 && proc_j == 0)
                continue;
            MPI_Recv(sub_C, block_dim*block_dim, c_type, proc_i*b_rows+proc_j, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (i = 0; i < block_dim; i++)
            {
                for (j = 0; j < block_dim; j++)
                {
                    C[(i + proc_i*block_dim)*dim + (j + proc_j*block_dim)] = sub_C[i*block_dim + j];
                }
            }
        }
        time_p = MPI_Wtime() - start;
        
        if (dim < 512 && compare_matrix(dim, Baseline, C) == -1)
        {
            printf("Matrix multiplication is wrong\n");
            if (dim <= 16)
            {
                printf("Matrix A: \n");
                print_matrix(dim, dim, A);
                printf("Matrix B: \n");
                print_matrix(dim, dim, B);
                printf("Matrix Baseline: \n");
                print_matrix(dim, dim, Baseline);
                printf("Matrix C: \n");
                print_matrix(dim, dim, C);
            }
            return -1;
        }
        MPI_Type_free(&row_vec);
        printf("For matrix size %d : Time for serial: %f ms \t Time for mpi_parallel: %f ms \n",dim, time_s*1000, time_p*1000);
    }
    else
    {
        sub_A = malloc((block_dim * dim) * sizeof(complex)); 
        sub_B = malloc((block_dim * dim) * sizeof(complex));
        sub_C = malloc((block_dim * block_dim) * sizeof(complex));

        MPI_Recv(sub_A, block_dim*dim, c_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(sub_B, block_dim*dim, c_type, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        multiply_naive(block_dim, dim, sub_A, dim, block_dim, sub_B, sub_C);

        MPI_Send(sub_C, block_dim*block_dim, c_type, 0, 9, MPI_COMM_WORLD);
    }
    MPI_Type_free(&c_type);
    MPI_Finalize();
    return 0;
}