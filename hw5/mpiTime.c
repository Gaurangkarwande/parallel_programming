#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_ITER 1000

void initialize_matrix(int dim, double mat[])
{
    int i;
    for ( i = 0; i < dim; i++)
        mat[i] = 3.14;
   return; 
}

int main(int argc, char *argv[])
{
    int dim = 2000;
    double start, time, *word, *recv_word;
    int rank, numprocs;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    word = malloc(dim * sizeof(double));
    initialize_matrix(dim, word);

    recv_word = malloc(dim * sizeof(double));
 
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        start = MPI_Wtime();

    for (int i = 0; i < NUM_ITER; i++)
    {
        if (rank == 0)
        {
            MPI_Send(word, dim, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
        if (rank == 1)
        {
            MPI_Recv(recv_word, dim, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0)
    {
        time = MPI_Wtime() - start;
        printf("The round-trip time for %d iterations and %d sized word was %f ms \n", NUM_ITER, dim, time*1000);
    }
    
    MPI_Finalize();
    return 0;
}