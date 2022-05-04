#include "main.h"

#define NUM_NODES 18772
#define NUM_EDGES 18772
#define SOURCE 9

int main(int argc, char *argv[])
{
    int num_nodes = NUM_NODES;
    int num_edges =  NUM_EDGES;
    int source = SOURCE;
    time_t start, time_s, time_p;
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

    // print_int_array(num_nodes, V);
    // print_int_array(num_nodes+1, I);
    // print_int_array(num_edges, E);
    // print_float_array(NUM_EDGES, W);

    start = time(NULL);
    bellman_sequential(source, num_nodes, num_edges, V, I, E, W, "./output/sequential_directed.txt");
    time_s = time(NULL) - start;

    read_int_array(V_file, V);
    read_int_array(I_file, I);
    read_int_array(E_file, E);
    read_float_array(W_file, W);

    start = time(NULL);
    bellman_parallel(source, num_nodes, num_edges, V, I, E, W, "./output/parallel_directed.txt");
    time_p = time(NULL) - start;

    printf("Time taken for serial: %ld ms, time taken for parallel: %ld ms \n", time_s*1000, time_p*1000);
    
    return 0;
}