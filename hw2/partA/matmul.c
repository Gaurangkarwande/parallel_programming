#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NUM_THREADS 8
#define BLOCK_SIZE 32
#define min(a,b) (((a)<(b))?(a):(b))

void initialize_matrix(int row, int col, double mat[], int is_result)
{
    srand(0);
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            if (is_result) {
                mat[i*col + j] = 0;
            } else {
                mat[i*col + j] = rand() % (row*col);
            } 
        }
    }
   return; 
}

int compare_matrix(int r1, int c1, double first[], int r2, int c2, double second[])
{
    for (int i = 0; i < r1; ++i)
        for (int j = 0; j < c2; ++j) 
            if (first[i*c2 + j] != second[i*c2 + j])
                return -1;
    return 0;
}

void multiply_serial_naive(int r1, int c1, double first[], int r2, int c2, double second[], double result[])
{
    int i,j,k;
    double first_ik;
    for (i = 0; i < r1; i++) {
        for (k = 0; k < c1; k++) {
            first_ik = first[i*c1 + k];
            for (j = 0; j < c2; j++) {
                result[i*c2 + j] += first_ik * second[k*c2 + j];
            }
        }
    }
}

void multiply_serial_block(int r1, int c1, double first[], int r2, int c2, double second[], double result[])
{
    int ii, jj, kk, i, j, k, first_ik, blk_size = BLOCK_SIZE;

	for(ii = 0; ii < r1; ii += blk_size) {
        for(kk = 0; kk < c1; kk += blk_size) {
            for(jj = 0; jj < c2; jj += blk_size) {
				for(i = ii; i < min(r1, ii+blk_size); i++) {
                    for(k = kk; k < min(c1, kk+blk_size); k++) {
                        first_ik = first[i*c1 + k];
                        for(j = jj; j < min(c2, jj+blk_size); j++) {
							result[i*c2 + j] += first_ik * second[k*c2 + j]; 
                        }
                    }
                }
            }
        }
    }											
}

void multiply_omp(int r1, int c1, double first[], int r2, int c2, double second[], double result[])
{
    int i,j,k;
    double first_ik;
    #pragma omp parallel for collapse(2) private(i,j,k) shared(first, second, result, r1, r2, c1, c2) num_threads(64)
        for (i = 0; i < r1; i++) {
            for (k = 0; k < c1; k++) {
                first_ik = first[i*c1 + k];
                for (j = 0; j < c2; j++) {
                    result[i*c2 + j] += first_ik * second[k*c2 + j];
                }
            }
        }
}

void multiply_omp_block(int r1, int c1, double first[], int r2, int c2, double second[], double result[])
{
    int ii, jj, kk, i, j, k, first_ik, blk_size = BLOCK_SIZE;

    #pragma omp parallel for private(i,j,k, ii, jj, kk) shared(first, second, result, r1, r2, c1, c2, blk_size) num_threads(32)
	for(ii = 0; ii < r1; ii += blk_size) {
        for(kk = 0; kk < c1; kk += blk_size) {
            for(jj = 0; jj < c2; jj += blk_size) {
				for(i = ii; i < min(r1, ii+blk_size); i++) {
                    for(k = kk; k < min(c1, kk+blk_size); k++) {
                        first_ik = first[i*c1 + k];
                        for(j = jj; j < min(c2, jj+blk_size); j++) {
							result[i*c2 + j] += first_ik * second[k*c2 + j]; 
                        }
                    }
                }
            }
        }
    }											
}

void multiply_omp_sections(int r1, int c1, double first[], int r2, int c2, double second[], double result[])
{
    int i,j,k;
    double first_ik = 0.0;
        for (i = 0; i < r1; i++) {
            #pragma omp parallel sections shared(first, second, result, r1, r2, c1, c2)
            {
            #pragma omp section
                for (k = 0; k < c1; k++) {
                    first_ik = first[i*c1 + k];
                    for (j = 0; j < c2; j++) {
                        result[i*c2 + j] += first_ik * second[k*c2 + j];
                    }
                }
            }
        }
}

int main()
{
    int dim = 1024;
    double start, time_s, time_p;
    double* A = malloc((dim * dim) * sizeof(double));
    double* B = malloc((dim * dim) * sizeof(double));
    double* C = malloc((dim * dim) * sizeof(double));
    double* Baseline = malloc((dim * dim) * sizeof(double));

    initialize_matrix(dim, dim, A, 0);
    initialize_matrix(dim, dim, B, 0);
    initialize_matrix(dim, dim, C, 1);
    initialize_matrix(dim, dim, Baseline, 1);

    start = omp_get_wtime();
    multiply_serial_naive(dim, dim, A, dim, dim, B, Baseline);
    time_s = omp_get_wtime() - start;

    start = omp_get_wtime();
    multiply_omp_sections(dim, dim, A, dim, dim, B, C);
    time_p = omp_get_wtime() - start;


    if (compare_matrix(dim, dim, Baseline, dim, dim, C) == -1)
    {
        printf("Matrix multiplication is wrong\n");
        return -1;
    }

    printf("Time for serial: %f ms \t Time for parallel: %f ms \n", time_s*1000, time_p*1000);
    return 0;
}