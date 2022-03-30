#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int r, i;
} complex;

void initialize_matrix(int dim, int mat[], int result)
{
    int i, j;
    for ( i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            mat[i*dim + j] = 0;
            if (i == j && !result)
                mat[i*dim + j] = 1;
        }
    }
   return; 
}

void transpose_matrix(int dim, int mat[], int result[])
{
    int i, j;
    for ( i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            result[j*dim + i] = mat[i*dim + j];
   return; 
}

void print_matrix(int row, int col, int mat[])
{
    int i, j;
    printf("\n");
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d \t", mat[i*col + j]);
        }
        printf("\n");
    }
    printf("\n");
   return; 
}

int compare_matrix(int dim, int first[], int second[])
{
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) 
            if (first[i*dim + j] != second[i*dim + j])
                return -1;
    return 0;
}

void multiply_naive(int first[], int second[], int result[], int dim)
{
    int i,j,k, sum;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            sum = 0;
            for (k = 0; k < dim; k++) {
                sum += first[i*dim + k] * second[k*dim + j];
            }
            result[i*dim + j] = sum;
        }
    }
}

void multiply_transpose(int first[], int second[], int result[], int block_dim, int dim)
{
    //Second matrix is transposed

    int i,j,k, sum;
    for (i = 0; i < block_dim; i++) {
        for (j = 0; j < block_dim; j++) {
            sum = 0;
            for (k = 0; k < dim; k++) {
                sum += first[i*dim + k] * second[j*dim + k];
            }
            result[i*block_dim + j] = sum;
        }
    }
}

void multiply_submatrix(int first[], int second[], int result[], int block_dim, int dim)
{
    //Second matrix is transposed

    int i,j,k, sum;
    for (i = 0; i < block_dim; i++) {
        for (j = 0; j < block_dim; j++) {
            sum = 0;
            for (k = 0; k < dim; k++) {
                sum += first[i*dim + k] * second[j*dim + k];
            }
            result[i*dim + j] = sum;
        }
    }
}

void multiply_mpi(int A[], int Bt[], int C[], int dim, int argc, char *argv[])
{
    int *sub_A, *sub_B, *sub_C;
    int rank, numprocs, block_dim, b_rows, b_cols, i, j, proc_i, proc_j;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    block_dim = dim/sqrt(numprocs);
    b_rows = b_cols = (int) sqrt(numprocs);

    if (rank == 0)  
    {
        sub_C = malloc((block_dim*block_dim) * sizeof(int));
        // print_matrix(dim, dim, C);
        for (i = 0; i < b_rows; i++)
        {
            for (j = 0; j < b_cols; j++)
            {
                if (i == 0 && j == 0)
                    continue;
                MPI_Send(A+i*block_dim*dim, block_dim*dim, MPI_INT, i*b_rows+j, 0, MPI_COMM_WORLD); //proc_id = i*b_rows+j
                MPI_Send(Bt+j*block_dim*dim, block_dim*dim, MPI_INT, i*b_rows+j, 1, MPI_COMM_WORLD);
                
            }
            
        }

        multiply_submatrix(A, Bt, C, block_dim, dim);
        // print_matrix(dim, dim, C);
        for (proc_i = 0; proc_i < b_rows; proc_i++)
        {
            for (proc_j = 0; proc_j < b_cols; proc_j++)
            {
                if (proc_i == 0 && proc_j == 0)
                    continue;
                MPI_Recv(sub_C, block_dim*block_dim, MPI_INT, proc_i*b_rows+proc_j, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (i = 0; i < block_dim; i++)
                {
                    for (j = 0; j < block_dim; j++)
                    {
                        C[(i + proc_i*block_dim)*dim + (j + proc_j*block_dim)] = sub_C[i*block_dim + j];
                    }
                }
            }
        }
        // print_matrix(dim, dim, C);
        
        
        
    }
    else
    {
        sub_A = malloc((block_dim * dim) * sizeof(int)); 
        sub_B = malloc((block_dim * dim) * sizeof(int));
        sub_C = malloc((block_dim * block_dim) * sizeof(int));

        MPI_Recv(sub_A, block_dim*dim, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(sub_B, block_dim*dim, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        multiply_transpose(sub_A, sub_B, sub_C, block_dim, dim);
        // print_matrix(block_dim, block_dim, sub_C);

        MPI_Send(sub_C, block_dim*block_dim, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}

int main(int argc, char *argv[])
{
    int dim = 1260;
    double start, time_s, time_p;
    int* A = malloc((dim * dim) * sizeof(int));
    int* B = malloc((dim * dim) * sizeof(int));
    int* C = malloc((dim * dim) * sizeof(int));
    int* Baseline = malloc((dim * dim) * sizeof(int));
    int* Bt = malloc((dim * dim) * sizeof(int));

    initialize_matrix(dim, A, 0);
    initialize_matrix(dim, B, 0);
    initialize_matrix(dim, C, 1 );
    initialize_matrix(dim, Baseline, 1);
    transpose_matrix(dim, B, Bt);
    
    start = MPI_Wtime();
    multiply_naive(A, B, Baseline, dim);
    time_s = MPI_Wtime() - start;

    int *sub_A, *sub_B, *sub_C;
    int rank, numprocs, block_dim, b_rows, b_cols, i, j, proc_i, proc_j;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    block_dim = dim/sqrt(numprocs);
    b_rows = b_cols = (int) sqrt(numprocs);

    if (rank == 0)  
    {
        start = MPI_Wtime();
        sub_C = malloc((block_dim*block_dim) * sizeof(int));
        // print_matrix(dim, dim, C);
        for (i = 0; i < b_rows; i++)
        {
            for (j = 0; j < b_cols; j++)
            {
                if (i == 0 && j == 0)
                    continue;
                MPI_Send(A+i*block_dim*dim, block_dim*dim, MPI_INT, i*b_rows+j, 0, MPI_COMM_WORLD);
                MPI_Send(Bt+j*block_dim*dim, block_dim*dim, MPI_INT, i*b_rows+j, 1, MPI_COMM_WORLD);
                
            }
            
        }

        multiply_submatrix(A, Bt, C, block_dim, dim);
        // print_matrix(dim, dim, C);
        for (proc_i = 0; proc_i < b_rows; proc_i++)
        {
            for (proc_j = 0; proc_j < b_cols; proc_j++)
            {
                if (proc_i == 0 && proc_j == 0)
                    continue;
                MPI_Recv(sub_C, block_dim*block_dim, MPI_INT, proc_i*b_rows+proc_j, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (i = 0; i < block_dim; i++)
                {
                    for (j = 0; j < block_dim; j++)
                    {
                        C[(i + proc_i*block_dim)*dim + (j + proc_j*block_dim)] = sub_C[i*block_dim + j];
                    }
                }
            }
        }
        time_p = MPI_Wtime() - start;
        
        if (compare_matrix(dim, Baseline, C) == -1)
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
        printf("Time for serial: %f ms \t Time for mpi_parallel: %f ms \n", time_s*1000, time_p*1000);
    }
    else
    {
        sub_A = malloc((block_dim * dim) * sizeof(int)); 
        sub_B = malloc((block_dim * dim) * sizeof(int));
        sub_C = malloc((block_dim * block_dim) * sizeof(int));

        MPI_Recv(sub_A, block_dim*dim, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(sub_B, block_dim*dim, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        multiply_transpose(sub_A, sub_B, sub_C, block_dim, dim);
        // print_matrix(block_dim, block_dim, sub_C);

        MPI_Send(sub_C, block_dim*block_dim, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}