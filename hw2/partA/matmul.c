#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NUM_THREADS 256
#define BLOCK_SIZE 32
#define min(a,b) (((a)<(b))?(a):(b))

void initialize_matrix(int row, int col, int mat[])
{
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            mat[i*col + j] = 0;
            if (i == j)
                mat[i*col + j] = 1;
        }
    }
   return; 
}

void print_matrix(int row, int col, int mat[])
{
    int i, j;
    for ( i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d \t", mat[i*col + j]);
        }
        printf("\n");
    }
   return; 
}

int compare_matrix(int r1, int c1, int first[], int r2, int c2, int second[])
{
    for (int i = 0; i < r1; ++i)
        for (int j = 0; j < c2; ++j) 
            if (first[i*c2 + j] != second[i*c2 + j])
                return -1;
    return 0;
}

void multiply_naive(int first[], int second[], int result[], int dim)
{
    int i,j,k, sum;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            sum = 0;
            for (k = 0; k < dim; k++) {
                sum += first[i*dim + k] * second[k*dim + j];
            }
            result[i*dim + j] = sum;
        }
    }
}

void multiply_serial_block(int first[], int second[], int result[], int dim)
{
    int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;

	for(ii = 0; ii < dim; ii += blk_size) {
        for(jj = 0; jj < dim; jj += blk_size) {
            for(kk = 0; kk < dim; kk += blk_size) {
				for(i = ii; i < min(dim, ii+blk_size); i++) {
                    for(j = jj; j < min(dim, jj+blk_size); j++) {
                        sum = 0;
                        for(k = kk; k < min(dim, kk+blk_size); k++) {
							sum += first[i*dim + k] * second[k*dim + j]; 
                        }
                        sum = result[i*dim + j];
                    }
                }
            }
        }
    }											
}

void multiply_omp(int first[], int second[], int result[], int dim)
{
    int i,j,k, sum;
    #pragma omp parallel for collapse(2) private(i,j,k, sum) shared(first, second, result)
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                sum = 0;
                for (k = 0; k < dim; k++) {
                    sum += first[i*dim + k] * second[k*dim + j];
                }
                result[i*dim + j] = sum;
            }
    }
}

void multiply_omp_block(int first[], int second[], int result[], int dim)
{
    int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;

    #pragma omp parallel for collapse(3) private(i,j,k, ii, jj, kk, sum) shared(first, second, result, blk_size)
        for(ii = 0; ii < dim; ii += blk_size) {
            for(jj = 0; jj < dim; jj += blk_size) {
                for(kk = 0; kk < dim; kk += blk_size) {
                    for(i = ii; i < min(dim, ii+blk_size); i++) {
                        for(j = jj; j < min(dim, jj+blk_size); j++) {
                            sum = 0;
                            for(k = kk; k < min(dim, kk+blk_size); k++) {
                                sum += first[i*dim + k] * second[k*dim + j]; 
                            }
                            sum = result[i*dim + j];
                        }
                    }
                }
            }
        }											
}

void multiply_omp_sections_block(int first[], int second[], int result[], int dim)
{
    int num_sections = 8;
    int section_size = BLOCK_SIZE/num_sections;
    #pragma omp parallel sections shared(first, second, result)
    {
        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*0; ii < section_size*(0+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*1; ii < section_size*(1+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*2; ii < section_size*(2+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*3; ii < section_size*(3+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*4; ii < section_size*(4+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*5; ii < section_size*(5+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*6; ii < section_size*(6+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }

        #pragma omp section
        {
            int ii, jj, kk, i, j, k,sum, blk_size = BLOCK_SIZE;
            for(ii = section_size*7; ii < section_size*(7+1); ii += blk_size) {
                for(jj = 0; jj < dim; jj += blk_size) {
                    for(kk = 0; kk < dim; kk += blk_size) {
                        for(i = ii; i < min(dim, ii+blk_size); i++) {
                            for(j = jj; j < min(dim, jj+blk_size); j++) {
                                sum = 0;
                                for(k = kk; k < min(dim, kk+blk_size); k++) {
                                    sum += first[i*dim + k] * second[k*dim + j]; 
                                }
                                sum = result[i*dim + j];
                            }
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel sections shared(first, second, result)
    {
        #pragma omp section
        {
            int i,j,k, sum;
            for (i = 0; i < section_size; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size; i < section_size*2; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }
    }
}

void multiply_omp_sections(int first[], int second[], int result[], int dim)
{
    int num_sections = 8;
    int section_size = dim/num_sections;
    #pragma omp parallel sections shared(first, second, result)
    {
        #pragma omp section
        {
            int i,j,k, sum;
            for (i = 0; i < section_size; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size; i < section_size*2; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*2; i < section_size*3; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*3; i < section_size*4; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*4; i < section_size*5; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*5; i < section_size*6; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*6; i < section_size*7; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

        #pragma omp section
        {
            int i,j,k, sum;
            for (i = section_size*7; i < section_size*8; i++) {
                for (j = 0; j < dim; j++) {
                    sum = 0;
                    for (k = 0; k < dim; k++) {
                        sum += first[i*dim + k] * second[k*dim + j];
                    }
                    result[i*dim + j] = sum;
                }
            }
        }

    }
}

int main()
{
    int dim = 2048;
    double start, time_s, time_p;
    int* A = malloc((dim * dim) * sizeof(int));
    int* B = malloc((dim * dim) * sizeof(int));
    int* C = malloc((dim * dim) * sizeof(int));
    int* Baseline = malloc((dim * dim) * sizeof(int));

    initialize_matrix(dim, dim, A);
    initialize_matrix(dim, dim, B);
    initialize_matrix(dim, dim, C);
    initialize_matrix(dim, dim, Baseline);

    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);

    // start = omp_get_wtime();
    // multiply_serial_block(A, B, Baseline, dim);
    // time_s = omp_get_wtime() - start;

    start = omp_get_wtime();
    multiply_omp_sections_block(A, B, C, dim);
    time_p = omp_get_wtime() - start;


    // if (compare_matrix(dim, dim, Baseline, dim, dim, C) == -1)
    // {
    //     printf("Matrix multiplication is wrong\n");
    //     if (dim <= 16)
    //     {
    //         printf("Matrix A: \n");
    //         print_matrix(dim, dim, A);
    //         printf("Matrix B: \n");
    //         print_matrix(dim, dim, B);
    //         printf("Matrix Baseline: \n");
    //         print_matrix(dim, dim, Baseline);
    //         printf("Matrix C: \n");
    //         print_matrix(dim, dim, C);
    //     }
    //     return -1;
    // }

    // printf("Time for serial: %f ms \t Time for loop parallel: %f ms \n", time_s*1000, time_p*1000);
    printf("Time for loop parallel: %f ms with %d threads and block size: %d \n", time_p*1000, NUM_THREADS, BLOCK_SIZE);
    printf("Time for section parallel: %f ms with %d threads \n", time_s*1000, NUM_THREADS);
    return 0;
}