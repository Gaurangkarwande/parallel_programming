#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>
#define BLOCK_SIZE 32
#define NUM_THREADS 16

void multiply_strassen(int A[], int B[], int C[], int dim);

struct thread_data
{
   int* A[4];
   int* B[4];
   int* P[7];
   int block_size
};

void initialize_matrix(int row, int col, int mat[])
{
    srand(0);
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

void multiply_naive(int r1, int c1, int first[], int r2, int c2, int second[], int result[])
{
    int i,j,k, sum;
    for (i = 0; i < r1; i++) {
        for (j = 0; j < c2; j++) {
            sum = 0;
            for (k = 0; k < c1; k++) {
                sum += first[i*c1 + k] * second[k*c2 + j];
            }
            result[i*c2 + j] = sum;
        }
    }
}

void split(int *C11, int *C12, int *C21, int *C22, int *C, int block_size) {
    // dim of sub_matrices C11.. is block_size, but dim of original matrix is block_size*2
	for(int i = 0; i < block_size; i++) {
		for(int j = 0; j < block_size; j++) {
			C11[i * block_size + j] = C[i * 2 * block_size + j];    //C[i][j]
			C12[i * block_size + j] = C[i * 2 * block_size + j + block_size];   //C[i][j+block_size]
			C21[i * block_size + j] = C[(i + block_size) * 2 * block_size + j]; //C[i+block_size][j]
			C22[i * block_size + j] = C[(i + block_size) * 2 * block_size + j + block_size];    //C[i+block_size][j+block_size]
		}
	}
}

/*
void copy_matrix_block(int *mat_block, int block_size, int *mat, int row, int col)
{
	for (int i = 0; i < block_size; i++) 
		mat_block[i] = &mat[(row + i) * block_size + col];
}
*/

void merge(int *C11, int *C12, int *C21, int *C22, int *C, int block_size) {
	for(int i = 0; i < block_size; i++) {
		for(int j = 0; j < block_size; j++) {
			C[i * 2 * block_size + j] = C11[i * block_size + j];
			C[i * 2 * block_size + j + block_size] = C12[i * block_size + j];
			C[(i + block_size) *2 * block_size + j] = C21[i * block_size + j];
			C[(i + block_size) * 2 * block_size + j + block_size] = C22[i * block_size + j];
		}
	}
    // print_matrix(block_size*2, block_size*2, C);
}

void compute_C(int *C11, int *C12, int *C21, int *C22, int *P1, int *P2, int *P3, int *P4, int *P5, int *P6, int *P7, int size)
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            C11[i * size + j] = P5[i * size + j] + P4[i * size + j] - P2[i * size + j] + P6[i * size + j];
            C12[i * size + j] = P1[i * size + j] + P2[i * size + j];
            C21[i * size + j] = P3[i * size + j] + P4[i * size + j];
            C22[i * size + j] = P5[i * size + j] + P1[i * size + j] - P3[i * size + j] - P7[i * size + j];
        }
    }
}

void add(int *A, int *B, int *C, int dim) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			C[i * dim + j] = A[i * dim + j] + B[i * dim + j];
		}
	}
}

void sub(int *A, int *B, int *C, int dim) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			C[i * dim + j] = A[i * dim + j] - B[i * dim + j];
		}
	}
}

void multiply_strassen_serial(int A[], int B[], int C[], int dim)
{
    if (dim <= BLOCK_SIZE){
        multiply_naive(dim, dim, A, dim, dim, B, C);
    }
    else
    {
        int block_size = dim / 2;
        int *A11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A22 = (int *) malloc((block_size*block_size) * sizeof(int));

        int *B11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B22 = (int *) malloc((block_size*block_size) * sizeof(int));

        int *C11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C22 = (int *) malloc((block_size*block_size) * sizeof(int));
        

        int *P1 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P2 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P3 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P4 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P5 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P6 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P7 = (int *) malloc((block_size*block_size) * sizeof(int));

        int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *add_result2 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *sub_result = (int *) malloc((block_size*block_size) * sizeof(int));

        split(A11, A12, A21, A22, A, block_size);
        split(B11, B12, B21, B22, B, block_size);

        //P1 = A11 * (B12 - B22)
        sub(B12, B22, sub_result, block_size);
        multiply_strassen_serial(A11, sub_result, P1, block_size);

        //P2 = (A11 + A12) * B22
        add(A11, A12, add_result1, block_size);
        multiply_strassen_serial(add_result1, B22, P2, block_size);

        //P3 = (A21 + A22) * B11
        add(A21, A22, add_result1, block_size);
        multiply_strassen_serial(add_result1, B11, P3, block_size);

        //P4 = (A21 + A22) * B11
        sub(B21, B11, sub_result, block_size);
        multiply_strassen_serial(A22, sub_result, P4, block_size);

        //P5 = (A11 + A22) * (B11 + B22)
        add(A11, A22, add_result1, block_size);
        add(B11, B22, add_result2, block_size);
        multiply_strassen_serial(add_result1, add_result2, P5, block_size);

        //P6 = (A12 - A22) * (B21 + B22)
        sub(A12, A22, sub_result, block_size);
        add(B21, B22, add_result1, block_size);
        multiply_strassen_serial(sub_result, add_result1, P6, block_size);

        //P7 = (A11 - A21) * (B11 + B12)
        sub(A11, A21, sub_result,block_size);
        add(B11, B12, add_result1, block_size);
        multiply_strassen_serial(sub_result, add_result1, P7, block_size);

        //deallocate A, B submatrices and add, sub results
        free(A11); free(A12); free(A21); free(A22);
        free(B11); free(B12); free(B21); free(B22);
        free(add_result1); free(add_result2); free(sub_result);

        //Compute C matrix blocks
        compute_C(C11, C12, C21, C22, P1, P2, P3, P4, P5, P6, P7, block_size);

        //deallocate P submatrices
        free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);

        //merge blocks

        merge(C11, C12, C21, C22, C, block_size);
        free(C11); free(C12); free(C21); free(C22);

    }
}



void* compute_P1(void *input)// (int B12[], int B22[], int A11[], int P1[], int block_size)
{
    int block_size = ((struct thread_data*)input)->block_size;
    int *A11 = ((struct thread_data*)input)->A[0];
    // int *A12 = ((struct thread_data*)input)->A[1];
    // int *A21 = ((struct thread_data*)input)->A[2];
    // int *A22 = ((struct thread_data*)input)->A[3];
    
    // int *B11 = ((struct thread_data*)input)->B[0];
    int *B12 = ((struct thread_data*)input)->B[1];
    // int *B21 = ((struct thread_data*)input)->B[2];
    int *B22 = ((struct thread_data*)input)->B[3];

    int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];

    int *sub_result = (int *) malloc((block_size*block_size) * sizeof(int));
    sub(B12, B22, sub_result, block_size);
    multiply_strassen(A11, sub_result, P1, block_size);
}

void* compute_P2(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    int *A11 = ((struct thread_data*)input)->A[0];
    int *A12 = ((struct thread_data*)input)->A[1];
    // int *A21 = ((struct thread_data*)input)->A[2];
    // int *A22 = ((struct thread_data*)input)->A[3];
    
    // int *B11 = ((struct thread_data*)input)->B[0];
    // int *B12 = ((struct thread_data*)input)->B[1];
    // int *B21 = ((struct thread_data*)input)->B[2];
    int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];

    int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
    add(A11, A12, add_result1, block_size);
    multiply_strassen(add_result1, B22, P2, block_size);
}

void* compute_P3(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    // int *A11 = ((struct thread_data*)input)->A[0];
    // int *A12 = ((struct thread_data*)input)->A[1];
    int *A21 = ((struct thread_data*)input)->A[2];
    int *A22 = ((struct thread_data*)input)->A[3];
    
    int *B11 = ((struct thread_data*)input)->B[0];
    // int *B12 = ((struct thread_data*)input)->B[1];
    // int *B21 = ((struct thread_data*)input)->B[2];
    // int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];

    int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
    add(A21, A22, add_result1, block_size);
    multiply_strassen(add_result1, B11, P3, block_size);
}

void* compute_P4(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    // int *A11 = ((struct thread_data*)input)->A[0];
    // int *A12 = ((struct thread_data*)input)->A[1];
    // int *A21 = ((struct thread_data*)input)->A[2];
    int *A22 = ((struct thread_data*)input)->A[3];
    
    int *B11 = ((struct thread_data*)input)->B[0];
    // int *B12 = ((struct thread_data*)input)->B[1];
    int *B21 = ((struct thread_data*)input)->B[2];
    // int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];

    int *sub_result = (int *) malloc((block_size*block_size) * sizeof(int));
    sub(B21, B11, sub_result, block_size);
    multiply_strassen(A22, sub_result, P4, block_size);
}

void* compute_P5(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    int *A11 = ((struct thread_data*)input)->A[0];
    // int *A12 = ((struct thread_data*)input)->A[1];
    // int *A21 = ((struct thread_data*)input)->A[2];
    int *A22 = ((struct thread_data*)input)->A[3];
    
    int *B11 = ((struct thread_data*)input)->B[0];
    // int *B12 = ((struct thread_data*)input)->B[1];
    // int *B21 = ((struct thread_data*)input)->B[2];
    int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];

    int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
    int *add_result2 = (int *) malloc((block_size*block_size) * sizeof(int));
    add(A11, A22, add_result1, block_size);
    add(B11, B22, add_result2, block_size);
    multiply_strassen(add_result1, add_result2, P5, block_size);
}

void* compute_P6(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    // int *A11 = ((struct thread_data*)input)->A[0];
    int *A12 = ((struct thread_data*)input)->A[1];
    // int *A21 = ((struct thread_data*)input)->A[2];
    int *A22 = ((struct thread_data*)input)->A[3];
    
    // int *B11 = ((struct thread_data*)input)->B[0];
    // int *B12 = ((struct thread_data*)input)->B[1];
    int *B21 = ((struct thread_data*)input)->B[2];
    int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    int *P6 = ((struct thread_data*)input)->P[5];
    // int *P7 = ((struct thread_data*)input)->P[6];
    
    int *sub_result = (int *) malloc((block_size*block_size) * sizeof(int));
    int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
    sub(A12, A22, sub_result, block_size);
    add(B21, B22, add_result1, block_size);
    multiply_strassen(sub_result, add_result1, P6, block_size);
}

void* compute_P7(void *input)
{
    int block_size = ((struct thread_data*)input)->block_size;
    int *A11 = ((struct thread_data*)input)->A[0];
    // int *A12 = ((struct thread_data*)input)->A[1];
    int *A21 = ((struct thread_data*)input)->A[2];
    // int *A22 = ((struct thread_data*)input)->A[3];
    
    int *B11 = ((struct thread_data*)input)->B[0];
    int *B12 = ((struct thread_data*)input)->B[1];
    // int *B21 = ((struct thread_data*)input)->B[2];
    // int *B22 = ((struct thread_data*)input)->B[3];

    // int *P1 = ((struct thread_data*)input)->P[0];
    // int *P2 = ((struct thread_data*)input)->P[1];
    // int *P3 = ((struct thread_data*)input)->P[2];
    // int *P4 = ((struct thread_data*)input)->P[3];
    // int *P5 = ((struct thread_data*)input)->P[4];
    // int *P6 = ((struct thread_data*)input)->P[5];
    int *P7 = ((struct thread_data*)input)->P[6];

    int *sub_result = (int *) malloc((block_size*block_size) * sizeof(int));
    int *add_result1 = (int *) malloc((block_size*block_size) * sizeof(int));
    sub(A11, A21, sub_result,block_size);
    add(B11, B12, add_result1, block_size);
    multiply_strassen(sub_result, add_result1, P7, block_size);
}

void multiply_strassen(int A[], int B[], int C[], int dim)
{
    if (dim <= BLOCK_SIZE){
        multiply_naive(dim, dim, A, dim, dim, B, C);
    }
    else
    {
        int block_size = dim / 2;
        int *A11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *A22 = (int *) malloc((block_size*block_size) * sizeof(int));

        int *B11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *B22 = (int *) malloc((block_size*block_size) * sizeof(int));

        int *C11 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C12 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C21 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *C22 = (int *) malloc((block_size*block_size) * sizeof(int));
        

        int *P1 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P2 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P3 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P4 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P5 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P6 = (int *) malloc((block_size*block_size) * sizeof(int));
        int *P7 = (int *) malloc((block_size*block_size) * sizeof(int));

        struct thread_data *data = (struct thread_data *)malloc(sizeof(struct thread_data));

        data->A[0] = A11;
        data->A[1] = A12;
        data->A[2] = A21;
        data->A[3] = A22;

        data->B[0] = B11;
        data->B[1] = B12;
        data->B[2] = B21;
        data->B[3] = B22;
        
        data->P[0] = P1;
        data->P[1] = P2;
        data->P[2] = P3;
        data->P[3] = P4;
        data->P[4] = P5;
        data->P[5] = P6;
        data->P[6] = P7;

        data->block_size = block_size;


        split(A11, A12, A21, A22, A, block_size);
        split(B11, B12, B21, B22, B, block_size);
        pthread_t threads[7];
        
        //P1 = A11 * (B12 - B22)
        pthread_create(&threads[0], NULL, compute_P1, (void *)data);

        //P2 = (A11 + A12) * B22
        pthread_create(&threads[1], NULL, compute_P2, (void *)data);

        //P3 = (A21 + A22) * B11
        pthread_create(&threads[2], NULL, compute_P3, (void *)data);

        //P4 = (A21 + A22) * B11
        pthread_create(&threads[3], NULL, compute_P4, (void *)data);

        //P5 = (A11 + A22) * (B11 + B22)
        pthread_create(&threads[4], NULL, compute_P5, (void *)data);

        //P6 = (A12 - A22) * (B21 + B22)
        pthread_create(&threads[5], NULL, compute_P6, (void *)data);

        //P7 = (A11 - A21) * (B11 + B12)
        pthread_create(&threads[6], NULL, compute_P7, (void *)data);;

        for (int tn = 0; tn < 7; tn++)
            pthread_join(threads[tn], NULL);
        
        // printf("Calculating P%d", p);
        //deallocate A, B submatrices and add, sub results
        // printf("\n Level %d: Merging blocks of size %d, into size %d by thread %d of %d\n", level, block_size, dim, omp_get_thread_num(), omp_get_num_threads());
        free(A11); free(A12); free(A21); free(A22);
        free(B11); free(B12); free(B21); free(B22);
        // free(add_result1); free(add_result2); free(sub_result);

        //Compute C matrix blocks
        compute_C(C11, C12, C21, C22, P1, P2, P3, P4, P5, P6, P7, block_size);

        //deallocate P submatrices
        free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);
        free(data);

        //merge blocks

        merge(C11, C12, C21, C22, C, block_size);
        free(C11); free(C12); free(C21); free(C22);

    }
}

void multiply_parallel(int A[], int B[], int C[], int dim)
{
    multiply_strassen(A, B, C, dim);
}

int main()
{
    int dim = 1024;
    double start, time_s, time_p;
    int* A = (int *) malloc((dim * dim) * sizeof(int));
    int* B = (int *) malloc((dim * dim) * sizeof(int));
    int* C = (int *) malloc((dim * dim) * sizeof(int));
    int* Baseline = (int *) malloc((dim * dim) * sizeof(int));

    printf("\nFor matrix of size %d x %d: \n", dim, dim);

    initialize_matrix(dim, dim, A);
    initialize_matrix(dim, dim, B);
    initialize_matrix(dim, dim, C);
    initialize_matrix(dim, dim, Baseline);

    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);
    

    start = omp_get_wtime();
    // multiply_naive(dim, dim, A, dim, dim, B, Baseline);
    multiply_strassen_serial(A, B, C, dim);
    time_s = omp_get_wtime() - start;

    start = omp_get_wtime();
    multiply_parallel(A, B, C, dim);
    time_p = omp_get_wtime() - start;


    if (compare_matrix(dim, dim, Baseline, dim, dim, C) == -1)
    {
        printf("Matrix multiplication is wrong\n");
        if (dim <= 16)
        {
            printf("Matrix A: \n");
            print_matrix(dim, dim, A);
            printf("Matrix B: \n");
            print_matrix(dim, dim, B);
            printf("Matrix Baseline: \n");
            print_matrix(dim, dim, Baseline);
            printf("Matrix C: \n");
            print_matrix(dim, dim, C);
        }
        return -1;
    }

    printf("Time for serial: %f ms \t Time for strassen_parallel: %f ms \n", time_s*1000, time_p*1000);
    return 0;
}