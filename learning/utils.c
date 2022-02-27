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