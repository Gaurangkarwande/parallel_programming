#include "../include/ssp_cuda.cuh"

__global__ void hellofunc()
{
    printf("Hello World\n");
}

void hello()
{
    hellofunc<<<1,1>>>();
    cudaDeviceSynchronize();
}

void printCudaDevice(){
    int dev = 0;
    cudaSetDevice(dev);
    cudaDeviceProp devProps;
    if (cudaGetDeviceProperties(&devProps, dev) == 0)
    {
        printf("****** Using device %d ***********\n", dev);
        printf("%s; global mem: %dB; compute v%d.%d; clock: %d kHz\n",
               devProps.name, (int)devProps.totalGlobalMem,
               (int)devProps.major, (int)devProps.minor,
               (int)devProps.clockRate);
        printf("Number of multiprocessors on device : %d\n", devProps.multiProcessorCount);
        printf("Maximum size of each dimension of a grid : %ld\n", devProps.maxGridSize);
        printf("Maximum size of each dimension of a block : %ld\n", devProps.maxThreadsDim);
        printf("Maximum number of threads per block : %d\n", devProps.maxThreadsPerBlock);
        //printf("Maximum number of resident blocks per multiprocessor : %d\n", devProps.maxBlocksPerMultiProcessor );
        printf("Maximum resident threads per multiprocessor : %d\n", devProps.maxThreadsPerMultiProcessor);
        printf("Shared memory available per block in bytes : %zu \n", devProps.sharedMemPerBlock );
        printf("Shared memory available per multiprocessor in bytes : %zu \n", devProps.sharedMemPerMultiprocessor );
        printf("Warp size in threads : %d \n", devProps.warpSize );
        printf("****** End of device stats ***********\n");
    }
}

__device__ __forceinline__ float atomicMinFloat (float * addr, float value) {
        float old;
        old = (value >= 0) ? __int_as_float(atomicMin((int *)addr, __float_as_int(value))) :
             __uint_as_float(atomicMax((unsigned int *)addr, __float_as_uint(value)));

        return old;
}

__global__ void init_distances_cuda(const int n, float D[], const int source)
{
    int index = threadIdx.x + blockDim.x * blockIdx.x;

    if (index < n) 
    {
        D[index] = __FLT_MAX__;
        if(index == source-1)
            D[index] = 0.0;
    }
}

__global__ void updateIndexOfEdges_cuda(int num_nodes, int num_edges, int V[], int E[]) 
{
    int l,r;
    unsigned int index = threadIdx.x + (blockDim.x * blockIdx.x);

    // This does binary search on the V array to find the index of each node in the Edge array (E) and replace the same with index
    // Based on the iterative binary search from : https://www.geeksforgeeks.org/binary-search/
    if (index < num_edges) {
        l=0; r=num_nodes-1;
        while (l <= r) {
            int m = l + (r - l) / 2;
            // Check if x is present at mid
            if (V[m] == E[index]) {
                E[index] = m;
                break;
            }
            // If x greater, ignore left half
            if (V[m] < E[index]) {
                l = m + 1;
            } else {        // If x is smaller, ignore right half
                r = m - 1;
            }
        }
    }
}

__global__ void relax(int num_nodes, float MAX_VAL, int V[], int I[], int E[], float W[], float D[], float Di[]) 
{
    unsigned int index = threadIdx.x + (blockDim.x * blockIdx.x);

    if (index < num_nodes) 
    {
        for (int j = I[index]; j < I[index + 1]; j++) 
        {
            // int u = V[index];
            // int v = V[E[j]];
            float w = W[j];
            float du = D[index];
            float dv = D[E[j]];
            float newDist = du + w;
            if (du == MAX_VAL)
                newDist = MAX_VAL;
            // printf("Index = %d, w=%f, du =%.1e, dv=%.1e,  -- du + w = %.1e\n", index, w, du , dv, du + w);

            if (newDist < dv)
                atomicMinFloat(&Di[E[j]], newDist);
        }
    }
}


__global__ void update_distance(int num_nodes, float D[], float Di[]) 
{
    unsigned int index = threadIdx.x + (blockDim.x * blockIdx.x);
    if (index < num_nodes) 
    {
        if (D[index] > Di[index])
            D[index] = Di[index];
        Di[index] = D[index];

    }
}

void bellman_parallel(int source, int num_nodes, int num_edges, int V[], int I[], int E[], float W[], const char* filename)
{
    printf("\n Running BellmanFord Parallel for graph with %d nodes and %d edges and source node %d \n", num_nodes, num_edges, source);
    float *D = (float *) malloc(num_nodes * sizeof(float));
    int *device_V, *device_I, *device_E;
    float *device_W, *device_D, *device_Di;

    cudaMalloc(&device_V, num_nodes*sizeof(int));
    cudaMalloc(&device_I, (num_nodes+1)*sizeof(int));
    cudaMalloc(&device_E, num_edges*sizeof(int));
    cudaMalloc(&device_W, num_edges*sizeof(float));

    cudaMalloc(&device_D, num_nodes*sizeof(float));
    cudaMalloc(&device_Di, num_nodes*sizeof(float));

    cudaMemcpy(device_V, V, num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(device_I, I, (num_nodes+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(device_E, E, num_edges*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(device_W, W, num_edges*sizeof(float), cudaMemcpyHostToDevice);

    int NUM_THREADS = 1024;
    int NUM_BLOCKS = ceil(num_edges/ (float) NUM_THREADS);
    init_distances_cuda<<<ceil(num_nodes/ (float) NUM_THREADS), NUM_THREADS>>>(num_nodes, device_D, source);
    init_distances_cuda<<<ceil(num_nodes/ (float) NUM_THREADS), NUM_THREADS>>>(num_nodes, device_Di, source);
    updateIndexOfEdges_cuda<<<ceil(num_edges/ (float) NUM_THREADS), NUM_THREADS>>>(num_nodes, num_edges, device_V, device_E);
    int *newE = (int *) malloc(num_edges * sizeof(int));

    for (int round = 1; round < num_nodes; round++)
    {
        relax<<<ceil(num_nodes/ (float) NUM_THREADS), NUM_THREADS>>>(num_nodes, __FLT_MAX__, device_V, device_I, device_E, device_W, device_D, device_Di);
        update_distance<<<ceil(num_nodes/ (float) NUM_THREADS), NUM_THREADS>>>(num_nodes, device_D, device_Di);
    }
    cudaDeviceSynchronize();
    cudaMemcpy(D, device_D, num_nodes*sizeof(float), cudaMemcpyDeviceToHost);
    save_to_file(num_nodes, V, D, filename);

    cudaFree(device_V);
    cudaFree(device_I);
    cudaFree(device_E);
    cudaFree(device_W);
    cudaFree(device_D);
    cudaFree(device_Di);
}