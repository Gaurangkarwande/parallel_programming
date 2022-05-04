#ifndef SSP_CUDA
#define SSP_CUDA

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include "utils.h"
// #include "kernels.cuh"

void hello();
void printCudaDevice();
void bellman_parallel(int source, int num_nodes, int num_edges, int V[], int I[], int E[], float W[], const char* filename);

#endif