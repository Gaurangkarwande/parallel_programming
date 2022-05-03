#ifndef SSP_UTILS
#define SSP_UTILS

#include <stdio.h>
#include <cstdlib>
#include <time.h>

int read_int_array(const char* filename, int A[]);
int read_float_array(const char* filename, float A[]);
void print_int_array(int n, int A[]);
void print_float_array(int n, float A[]);
void init_distances(int n, float D[]);
void save_to_file(int n, int V[], float D[], const char* filename);

#endif