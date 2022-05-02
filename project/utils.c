#include "utils.h"

int read_int_array(const char* filename, int A[])
{
    int i = 0;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 1;
    while (fscanf(pf, "%d", &A[i]) == 1)
    {
        i++;
    }
    fclose (pf);
    return 0;
}

int read_float_array(const char* filename, float A[])
{
    int i = 0;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 1;
    while (fscanf(pf, "%f", &A[i]) == 1)
    {
        i++;
    }
    fclose (pf);
    return 0;
}

void print_int_array(int n, int A[])
{
    for (int i = 0; i < n; i++)
        printf("%d\t", A[i]);
    printf("\n");
}

void print_float_array(int n, float A[])
{
    for (int i = 0; i < n; i++)
        printf("%.1f\t", A[i]);
    printf("\n");
}

void init_distances(int n, float D[])
{
    for (int i = 0; i < n; i++)
    {
        D[i] = __FLT_MAX__;
    }
    
}