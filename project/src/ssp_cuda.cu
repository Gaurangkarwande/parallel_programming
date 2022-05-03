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

// int main()
// {
//     printCudaDevice();
//     return 0;
// }