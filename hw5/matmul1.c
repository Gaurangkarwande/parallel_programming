#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void allocate_matrix(int row, int col, double ***mat)
{
    double *p = (double *) malloc((row * col) * sizeof(double*));
    *mat = (double **) malloc(row * sizeof(double*));
    for (int i = 0; i < row; i++)
        (*mat)[i] = &(p[i * col]);      //parenthesis is important
}

void initialize_matrix(int dim, double** mat, int result)
{
    srand(0);
    int i, j;
    for ( i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            mat[i][j] = 0;
            if (i == j && !result)
                mat[i][j] = 1;
        }
    }
   return; 
}

void transpose_matrix(int dim, double** mat, double** result)
{
    int i, j;
    for ( i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            result[j][i] = mat[i][j];
   return; 
}

void print_matrix(int row, int col, double** mat)
{
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%f \t", mat[i][j]);
        }
        printf("\n");
    }
   return; 
}

int compare_matrix(int r1, int c1, double** first, int r2, int c2, double** second)
{
    for (int i = 0; i < r1; ++i)
        for (int j = 0; j < c2; ++j) 
            if (first[i][j] != second[i][j])
                return -1;
    return 0;
}

void multiply_naive(int r1, int c1, double** first, int r2, int c2, double** second, double** result)
{
    int i,j,k, sum;
    for (i = 0; i < r1; i++) {
        for (j = 0; j < c2; j++) {
            sum = 0;
            for (k = 0; k < c1; k++) {
                sum += first[i][k] * second[k][j];
            }
            result[i][j] += sum;
        }
    }
}

int main(int argc, char *argv[])
{
    int dim = 420;
    double start, time_s, time_p;
    double **A, **B, **C, **Baseline;
    allocate_matrix(dim, dim, &A);
    allocate_matrix(dim, dim, &B);
    allocate_matrix(dim, dim, &C);
    allocate_matrix(dim, dim, &Baseline);

    initialize_matrix(dim, A, 0);
    initialize_matrix(dim, B, 0);
    initialize_matrix(dim, C, 1 );
    initialize_matrix(dim, Baseline, 1);
    
    start = MPI_Wtime();
    multiply_naive(dim, dim, A, dim, dim, B, Baseline);
    time_s = MPI_Wtime() - start;

    double **sub_A, **sub_B, **sub_C;
    int rank, numprocs, block_dim, b_rows, b_cols, cart_rank, cart_coords[2], left, right, up, down, i, j, k;
    
    MPI_Comm cart_comm;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (rank == 0)
        start = MPI_Wtime();

    block_dim = dim/sqrt(numprocs);
    b_rows = b_cols = (int) sqrt(numprocs);

    int cart_dim[2] = {b_rows, b_cols};
    int cart_period[2] = {1,1};

    MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dim, cart_period, 1, &cart_comm);
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 2, cart_coords);

    allocate_matrix(block_dim, block_dim, &sub_A);
    allocate_matrix(block_dim, block_dim, &sub_B);
    allocate_matrix(block_dim, block_dim, &sub_C);
    initialize_matrix(block_dim, sub_C, 1);

    int globalSize[2] = { dim, dim };
	int localSize[2] = { block_dim, block_dim };
	int starts[2] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, block_dim * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

	int* sendCounts = (int*)malloc(sizeof(int) * numprocs);
	int* displacements = (int*)malloc(sizeof(int) * numprocs);

	if (rank == 0) {
		for (int i = 0; i < numprocs; i++) {
			sendCounts[i] = 1;
		}
		int disp = 0;
		for (i = 0; i < b_rows; i++) {
			for (j = 0; j < b_cols; j++) {
				displacements[i * b_cols + j] = disp;
				disp += 1;
			}
			disp += (block_dim - 1)* b_rows;
		}
	}

	MPI_Scatterv(&(A[0][0]), sendCounts, displacements, subarrtype, &(sub_A[0][0]), block_dim*block_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(B[0][0]), sendCounts, displacements, subarrtype, &(sub_B[0][0]), block_dim*block_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Cart_coords(cart_comm, rank, 2, cart_coords);

    MPI_Cart_shift(cart_comm, 1, cart_coords[0], &left, &right);
    MPI_Sendrecv_replace(&(sub_A[0][0]), block_dim*block_dim, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);

    MPI_Cart_shift(cart_comm, 0, cart_coords[1], &up, &down);
    MPI_Sendrecv_replace(&(sub_B[0][0]), block_dim*block_dim, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
    

    for (k = 0; k < b_rows; k++)
    {
        multiply_naive(block_dim, block_dim, sub_A, block_dim, block_dim, sub_B, sub_C);
        
        MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
        MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
        MPI_Sendrecv_replace(&(sub_A[0][0]), block_dim*block_dim, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(&(sub_B[0][0]), block_dim*block_dim, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(cart_comm);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gatherv(&(sub_C[0][0]), block_dim*block_dim, MPI_DOUBLE, &(C[0][0]), sendCounts, displacements, subarrtype, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        time_p = MPI_Wtime() - start;
        if (dim < 512 && compare_matrix(dim, dim, Baseline, dim, dim, C) == -1)        //setting dim < 512 condition for experimentation
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
        printf("For matrix size %d, num procs %d : Time for serial: %f ms \t Time for mpi_parallel: %f ms \n",dim, numprocs, time_s*1000, time_p*1000);
    }

    MPI_Type_free(&subarrtype);

    MPI_Finalize();
    return 0;
}