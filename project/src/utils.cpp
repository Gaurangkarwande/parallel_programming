#include "../include/ssp_sequential.h"

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
        printf("%.1e\t", A[i]);
    printf("\n");
}

void init_distances(int n, float D[])
{
    for (int i = 0; i < n; i++)
        D[i] = __FLT_MAX__;   
}

void save_to_file(int n, int V[], float D[], const char* filename)
{
    printf("Saving result to file %s\n", filename);
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "Node \t:\tDistance\n");
    for (int i = 0; i < n; i++)
        fprintf(fp, "%d \t:\t%.1e\n", V[i], D[i]);
    fclose(fp);
}
