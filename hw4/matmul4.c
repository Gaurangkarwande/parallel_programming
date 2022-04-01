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

int main(int argc, char *argv[])
{
    int dim = 1260;
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

    complex *sub_A, *sub_B, *sub_C, *mat_C;
    int rank, numprocs, block_dim, b_rows, b_cols, i, j, k, proc_i, proc_j, row_rank, col_rank, numprocs_row, numprocs_col;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (rank == 0)
        start = MPI_Wtime();

    MPI_Datatype c_type;
    MPI_Type_contiguous(2, MPI_INT, &c_type);
    MPI_Type_commit(&c_type);

    block_dim = dim/sqrt(numprocs);
    b_rows = b_cols = (int) sqrt(numprocs);

    MPI_Datatype row_vec;
    MPI_Type_vector(dim, 1, 1, c_type, &row_vec);
    MPI_Type_commit(&row_vec);

    MPI_Comm_split(MPI_COMM_WORLD, rank/b_rows, rank, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, rank%b_cols, rank, &col_comm);

    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &numprocs_row);

    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &numprocs_col);

    sub_A = malloc((block_dim * dim) * sizeof(complex)); 
    sub_B = malloc((block_dim * dim) * sizeof(complex));
    sub_C = malloc((block_dim * block_dim) * sizeof(complex));

    proc_i = (int) rank / b_rows;
    proc_j = rank % b_cols;


    MPI_Scatter(A, block_dim, row_vec, sub_A, block_dim, row_vec, 0, col_comm );

    MPI_Scatter(Bt, block_dim, row_vec, sub_B, block_dim, row_vec, 0, row_comm);

    // printf("%d\n", rank);
    // print_matrix(block_dim, dim, sub_A);

    multiply_transpose(sub_A, sub_B, sub_C, block_dim, dim);
    
    // MPI_Barrier(col_comm);
    // MPI_Barrier(row_comm);
    // MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        mat_C = malloc((dim * dim) * sizeof(complex));

    MPI_Gather(sub_C, block_dim*block_dim, c_type, mat_C, block_dim*block_dim, c_type, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (k = 0; k < numprocs; k++)
        {
            proc_i = (int) k/b_rows;
            proc_j = k % b_cols;
            for (i = 0; i < block_dim; i++)
            {
                for (j = 0; j < block_dim; j++)
                {
                    C[(i + proc_i*block_dim)*dim + (j + proc_j*block_dim)] = mat_C[k*block_dim*block_dim + i*block_dim + j];
                }
            }
        }
        time_p = MPI_Wtime() - start;
        if (dim < 512 && compare_matrix(dim, Baseline, C) == -1)        //setting dim < 512 condition for experimentation
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
        }
        printf("Time for serial: %f ms \t Time for mpi_parallel: %f ms \n", time_s*1000, time_p*1000);
    }
    MPI_Type_free(&row_vec);
    MPI_Type_free(&c_type);
    MPI_Finalize();
    return 0;
}