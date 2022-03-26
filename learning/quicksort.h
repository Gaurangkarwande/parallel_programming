#ifndef QUICKSORT_H
#define QUICKSORT_H


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "utils.h"

int partition(int low, int high, int arr[]);
int partition_r(int low, int high, int arr[]);
void quicksort_ser(int low, int high, int arr[]);


#endif

