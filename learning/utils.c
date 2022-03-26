#include "utils.h"

void initialize_1Darray(int dim, int arr[])
{
    srand(345);
    int i;
    for ( i = 0; i < dim; i++) {
        arr[i] = rand() % dim;
    }
   return; 
}

int uniform_distribution(int rangeLow, int rangeHigh) {
    double myRand = rand()/(1.0 + RAND_MAX); 
    int range = rangeHigh - rangeLow + 1;
    int myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

void initialize_1Darray_inRange(int dim, int arr[], int rangeLow, int rangeHigh)
{
    srand(345);
    int i;
    for ( i = 0; i < dim; i++) {
        arr[i] = uniform_distribution(rangeLow, rangeHigh);
    }
   return; 
}

void print_1Darray(int dim, int arr[])
{
    int i;
    printf("Array is: \n");
    for ( i = 0; i < dim; i++) {
        printf("%d, ", arr[i]);
    }
    printf("\n");
   return; 
}


void swap(int *a, int *b)
{
   int temp = *a;
   *a = *b;
   *b = temp;
}