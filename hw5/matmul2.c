#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define iDIM 0
#define jDIM 1
#define kDIM 2

struct mesh_info {
	MPI_Comm mesh3d;
	MPI_Comm mesh_ik, ring_j;
	MPI_Comm mesh_jk, ring_i;	
	MPI_Comm mesh_ij, ring_k;

	int myrank;
	int coords[3]; // my coordinates in the 3D mesh

	int numprocs; // the number of processors in our cluster
	
	int q;	// the number of processors in each dimension of the 3D-Mesh 
		// each processor will receive two blocks of (n/q)*(n/q) elements

	int n; // the size of the matrix dimension
    int block_dim; // the size of the block dimension

	int k; // the power
};

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

void create_topology( struct mesh_info *info , int dim) {

	// MPI Initialization
	MPI_Comm_size(MPI_COMM_WORLD, &info->numprocs);	

    info->n = dim;
	info->q = (int) cbrt(info->numprocs);
    info->block_dim = (int) info->n/info->q;

	// Create the Topology which we will be using.
	int *dims = malloc( sizeof(int) * 3 ); // 3 dimensions
	int *periods = malloc( sizeof(int) * 3); // wraparound
	dims[iDIM] = dims[jDIM] = dims[kDIM] = info->q;
	periods[iDIM] = periods[jDIM] = periods[kDIM] = 1;
	
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &info->mesh3d);

	MPI_Comm_rank( info->mesh3d, &info->myrank);
	MPI_Cart_coords( info->mesh3d, info->myrank, 3, info->coords);

	// the i-k plane
	dims[iDIM] = dims[kDIM] = 1;
	dims[jDIM] = 0;

	MPI_Cart_sub( info->mesh3d, dims, &info->mesh_ik );

	// the j ring
	dims[iDIM] = dims[kDIM] = 0;
	dims[jDIM] = 1;
	
	MPI_Cart_sub( info->mesh3d, dims, &info->ring_j);

	// the j-k plane
	dims[jDIM] = dims[kDIM] = 1;
	dims[iDIM] = 0;

	MPI_Cart_sub( info->mesh3d, dims, &info->mesh_jk);

	//the i ring
	dims[jDIM] = dims[kDIM] = 0;
	dims[iDIM] = 1;

	MPI_Cart_sub( info->mesh3d, dims, &info->ring_i);

	// the i-j plane 
	dims[iDIM] = dims[jDIM] = 1;
	dims[kDIM] = 0;
	
	MPI_Cart_sub( info->mesh3d, dims, &info->mesh_ij);

	// the k rings
	dims[iDIM] = dims[jDIM] = 0;
	dims[kDIM] = 1;
    
	MPI_Cart_sub( info->mesh3d, dims, &info->ring_k);
}

/*
   Distribute the left or right hand matrix.
*/
void distribute_matrix(struct mesh_info *info, double** Mat, double** sub_Mat, int sendCounts[], int displacements[], int ringdim, int bycolumn) 
{   
    int globalSize[2] = { info->n, info->n };
	int localSize[2] = { info->block_dim, info->block_dim };
	int starts[2] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, info->block_dim * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

	MPI_Comm mesh, ring;	
	if (ringdim == jDIM) {
		mesh = info->mesh_ik;
		ring = info->ring_j;
	} else {
		mesh = info->mesh_jk;
		ring = info->ring_i;
	}

	if( info->coords[ringdim] == 0 ) {
        MPI_Scatterv(&(Mat[0][0]), sendCounts, displacements, subarrtype, &(sub_Mat[0][0]), info->block_dim*info->block_dim, MPI_DOUBLE, 0, mesh);
	}

    MPI_Bcast(&(sub_Mat[0][0]), info->block_dim*info->block_dim, MPI_DOUBLE, 0, ring);
    MPI_Type_free(&subarrtype);
}

void gather_blocks(struct mesh_info *info, double** Mat, double** sub_Mat, int sendCounts[], int displacements[])
{
    int globalSize[2] = { info->n, info->n };
	int localSize[2] = { info->block_dim, info->block_dim };
	int starts[2] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, info->block_dim * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

    MPI_Gatherv(&(sub_Mat[0][0]), info->block_dim*info->block_dim, MPI_DOUBLE, &(Mat[0][0]), sendCounts, displacements, subarrtype, 0, info->mesh_ij);

    MPI_Type_free(&subarrtype);
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

    double **sub_A, **sub_B, **sub_C, **sub_C_reduce;
    int i, j;
    struct mesh_info info;

    MPI_Init(&argc, &argv);	
    create_topology(&info, dim);

    allocate_matrix(info.block_dim, info.block_dim, &sub_A);
    allocate_matrix(info.block_dim, info.block_dim, &sub_B);
    allocate_matrix(info.block_dim, info.block_dim, &sub_C);
    allocate_matrix(info.block_dim, info.block_dim, &sub_C_reduce);

    initialize_matrix(info.block_dim, sub_C, 1);
    initialize_matrix(info.block_dim, sub_C_reduce, 1);
    
    int* sendCounts = (int*)malloc(sizeof(int) * info.q*info.q);
	int* displacements = (int*)malloc(sizeof(int) * info.q*info.q);

	if (info.myrank == 0) {
		for (i = 0; i < info.q*info.q; i++) {
			sendCounts[i] = 1;
		}
		int disp = 0;
		for (i = 0; i < info.q; i++) {
			for (j = 0; j < info.q; j++) {
				displacements[i * info.q + j] = disp;
				disp += 1;
			}
			disp += (info.block_dim - 1)* info.q;
		}
	}
    if (info.myrank == 0)
        start = MPI_Wtime();

    distribute_matrix(&info, A, sub_A, sendCounts, displacements, jDIM, 0);
    distribute_matrix(&info, B, sub_B, sendCounts, displacements, iDIM, 0);

    multiply_naive(info.block_dim, info.block_dim, sub_A, info.block_dim, info.block_dim, sub_B, sub_C);

    MPI_Reduce( (&sub_C[0][0]), &(sub_C_reduce[0][0]), info.block_dim*info.block_dim, MPI_DOUBLE, MPI_SUM, 0, info.ring_k);

    MPI_Barrier(MPI_COMM_WORLD);
    
    if (info.coords[kDIM] == 0)
    {
        gather_blocks(&info, C, sub_C_reduce, sendCounts, displacements);
    }
    
    if (info.myrank == 0)
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
        printf("For matrix size %d, num procs %d : Time for serial: %f ms \t Time for mpi_parallel: %f ms \n",dim, info.numprocs, time_s*1000, time_p*1000);
    }

    MPI_Finalize();
    return 0;
}