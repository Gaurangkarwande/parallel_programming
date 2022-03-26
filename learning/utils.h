#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

void initialize_1Darray(int dim, int arr[]);
void print_1Darray(int dim, int arr[]);
void swap(int *a, int *b);
void initialize_1Darray_inRange(int dim, int arr[], int rangeLow, int rangeHigh);

#endif