#include "main.h"

#define NUM_NODES 5
#define NUM_EDGES 9
#define SOURCE 5

int main(int argc, char *argv[])
{
    int num_nodes = NUM_NODES;
    int num_edges =  NUM_EDGES;
    int source = SOURCE;
    int *V = (int*) malloc(num_nodes * sizeof(int));
    int *I = (int*) malloc((num_nodes+1) * sizeof(int));
    int *E = (int*) malloc(num_edges * sizeof(int));
    float *W = (float*) malloc(num_edges * sizeof(float));
    
    char V_file[] = "./data/sample_V.txt";
    char I_file[] = "./data/sample_I.txt";
    char E_file[] = "./data/sample_E.txt";
    char W_file[] = "./data/sample_W.txt";

    read_int_array(V_file, V);
    read_int_array(I_file, I);
    read_int_array(E_file, E);
    read_float_array(W_file, W);

    // print_int_array(num_nodes, V);
    // print_int_array(num_nodes+1, I);
    // print_int_array(num_edges, E);
    // print_float_array(NUM_EDGES, W);

    bellman_sequential(source, num_nodes, num_edges, V, I, E, W, "./output/sample_directed.txt");
    printCudaDevice();
    
    return 0;
}