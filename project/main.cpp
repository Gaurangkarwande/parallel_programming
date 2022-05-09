#include "main.h"

#define NUM_NODES 18772
#define NUM_EDGES 18772
#define SOURCE 4

int main(int argc, char *argv[])
{
    int num_nodes = NUM_NODES;
    int num_edges =  NUM_EDGES;
    int source = SOURCE;
    struct timeval start, end;
    double time_s, time_p;
    int *V = (int*) malloc(num_nodes * sizeof(int));
    int *I = (int*) malloc((num_nodes+1) * sizeof(int));
    int *E = (int*) malloc(num_edges * sizeof(int));
    float *W = (float*) malloc(num_edges * sizeof(float));
    
    char V_file[] = "./data/Graph_directed_V.txt";
    char I_file[] = "./data/Graph_directed_I.txt";
    char E_file[] = "./data/Graph_directed_E.txt";
    char W_file[] = "./data/Graph_directed_W.txt";

    read_int_array(V_file, V);
    read_int_array(I_file, I);
    read_int_array(E_file, E);
    read_float_array(W_file, W);


    gettimeofday(&start, NULL);
    bellman_sequential(source, num_nodes, num_edges, V, I, E, W, "./output/sequential_directed.txt");
    gettimeofday(&end, NULL);
    time_s = (end.tv_sec - start.tv_sec) * 1e6;
    time_s = (time_s + (end.tv_usec - start.tv_usec)) * 1e-3;

    read_int_array(V_file, V);
    read_int_array(I_file, I);
    read_int_array(E_file, E);
    read_float_array(W_file, W);

    gettimeofday(&start, NULL);
    bellman_parallel(source, num_nodes, num_edges, V, I, E, W, "./output/parallel_directed.txt");
    gettimeofday(&end, NULL);
    time_p = (end.tv_sec - start.tv_sec) * 1e6;
    time_p = (time_p + (end.tv_usec - start.tv_usec)) * 1e-3;

    printf("Time taken for serial: %f ms, time taken for parallel: %f ms \n", time_s, time_p);
    
    return 0;
}