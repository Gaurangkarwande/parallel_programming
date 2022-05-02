#include "ssp_sequential.h"


void updateIndexOfEdges(int num_nodes, int num_edges, int V[], int E[])
{
    int l,r;
    for (int index = 0; index < num_edges; index++) {
        //cout << "Updating index of  E[index] " <<  E[index] << endl;
        l=0; r=num_nodes-1;
        while (l <= r) {
            int m = l + (r - l) / 2;
            // Check if x is present at mid
            if (V[m] == E[index]) {
                E[index] = m;
                break;
            }
            // If x greater, ignore left half
            if (V[m] < E[index]) {
                l = m + 1;
            } else {        // If x is smaller, ignore right half
                r = m - 1;
            }
        }
        //cout << "index of  E[index] " <<  E[index] << endl;
    }
}

void bellman_sequential(int source, int num_nodes, int num_edges, int V[], int I[], int E[], float W[])
{
    printf("\n Running BellmanFord Sequential for graph with %d nodes and %d edges and source node %d \n", num_nodes, num_edges, source);
    float *D = malloc(num_nodes * sizeof(float));
    int *pred = malloc(num_nodes * sizeof(int));
    

    init_distances(num_nodes, D);
    updateIndexOfEdges(num_nodes, num_edges, V, E);

    D[source-1] = 0;
    pred[source-1] = 0;


    // Bellman ford
    for (int round = 1; round < num_nodes; round++) {
        for (int i = 0; i < num_nodes ; i++) {
            for (int j = I[i]; j < I[i + 1]; j++) {
                int u = V[i];
                int v = V[E[j]];
                float w = W[j];
                float du = D[i];
                float dv = D[E[j]];
                if (du + w < dv) {
                    //cout<< "Relaxing edge (" << u << ", " << v << ") current dist =" << dv << ", new dist =" << du + w << endl;
                    D[E[j]] = du + w;
                    pred[E[j]] = u;
                }
            }
        }
    }

    print_float_array(num_nodes, D);
    
}