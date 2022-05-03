#ifndef SSP_SEQUENTIAL
#define SSP_SEQUENTIAL

#include <stdio.h>
#include <cstdlib>
#include "utils.h"

void bellman_sequential(int source, int num_nodes, int num_edges, int V[], int I[], int E[], float W[], const char* filename);

#endif