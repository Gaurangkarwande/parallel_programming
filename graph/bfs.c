#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>


typedef struct {
    int val, level;
} node;

int readmatrix(size_t rows, size_t cols, int mat[], const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 1;

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%d", &mat[i*cols + j]);
    }


    fclose (pf);
    return 0;
}
void print_matrix(int row, int col, int mat[])
{
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d \t", mat[i*col + j]);
        }
        printf("\n");
    }
	printf("\n");
   return; 
}

void initialize_nodes(int num_nodes, node graph[])
{
	for(int i = 0; i < num_nodes; i++)
	{
		graph[i].val = 0;
		graph[i].level = -1;
	}
}

void transpose_matrix(int dim, int mat[], int result[])
{
    int i, j;
    for ( i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            result[j*dim + i] = mat[i*dim + j];
   return; 
}

void multiply_matvec(int row, int col, int mat[], int vec[], int result[], int rank, int visited[])
{
	int sum;
    int i,j;
    for (i = 0; i < row; i++) {
		if(visited[row*rank + i] == 1)		//row = nodes_per_proc
			continue;
		sum = 0;
        for (j = 0; j < col; j++) {
            sum += mat[i*col + j] * vec[j];
        }
		sum = sum > 0 ? 1: 0;
		result[i] = sum;
    }
}

int main(int argc, char *argv[])
{
	//Variables and Initializations
	int num_nodes = 4;
	node *graph;
	int *A, *At, *sub_At, *visited, *frontier, *sub_next_frontier;

	graph = malloc(num_nodes * sizeof(node));
	A = malloc((num_nodes * num_nodes) * sizeof(int));
	At = malloc((num_nodes * num_nodes) * sizeof(int));
	visited = malloc(num_nodes * sizeof(int));
	frontier = malloc(num_nodes * sizeof(int));

	//MPI Code
	int num_procs, rank, i, nodes_per_proc, k;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	nodes_per_proc = num_nodes/num_procs;
	sub_At = malloc((nodes_per_proc * num_nodes) * sizeof(int));
	sub_next_frontier = malloc((nodes_per_proc) * sizeof(int));

	if(rank == 0)
	{
        if(readmatrix(num_nodes, num_nodes, A, "mat.dat") == 1)
        {
            printf("File not found \n");
            MPI_Finalize();
            return 1;
        }
		transpose_matrix(num_nodes, A, At);
		printf("The transposed adjacency matrix is: \n");
        print_matrix(num_nodes, num_nodes, At);

		for(int i = 0; i < num_nodes; i++)
		{
			visited[i] = 0;
			frontier[i] = 0;
		}

		initialize_nodes(num_nodes, graph);

		visited[0] = 1;
		graph[0].val = 1;
		graph[0].level = 0;
		frontier[0] = 1;

		printf("The initial frontier is: \n");
        print_matrix(num_nodes, 1, frontier);
	}
	//setup
	MPI_Scatter(At, nodes_per_proc * num_nodes, MPI_INT, sub_At, nodes_per_proc * num_nodes, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(frontier, num_nodes, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (i = 0; i < num_nodes; i++)
	{
		for (k = 0; k < nodes_per_proc; k++)
		{
			sub_next_frontier[k] = 0;
		}
		
		//compute
		multiply_matvec(nodes_per_proc, num_nodes, sub_At, frontier, sub_next_frontier, rank, visited);

		// now each proc has the current/next nodes to work upon. Can do processing of each node as per the problem here

		//communicate
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allgather(sub_next_frontier, nodes_per_proc, MPI_INT, frontier, nodes_per_proc, MPI_INT, MPI_COMM_WORLD);
		if (rank == 0)
		{
			for (k = 0; k < num_nodes; k++)
			{
				if (frontier[k] == 1)
				{
					visited[k] = 1;
					graph[k].val++;
					graph[k].level = i + 1;
				}
			}
		}
		//communicate
		MPI_Bcast(visited, num_nodes, MPI_INT, 0, MPI_COMM_WORLD);
	}

	if(rank == 0)
	{
        printf("\n BFS Traversal from node: 0 \n ");
		for(int i = 0; i < num_nodes ; i++)
		{
			printf("Node: %d, Value: %d, Level: %d \n", i, graph[i].val, graph[i].level);
		}
	}


	MPI_Finalize();
	return 0;
}