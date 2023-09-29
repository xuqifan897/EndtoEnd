#include <iostream>
#include <cuda_runtime.h>

int main() {
    // Initialize CUDA
    cudaFree(0);

    // Get device properties
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0); // Assumes device 0

    // Print resource limits
    std::cout << "Device Name: " << deviceProp.name << std::endl;
    std::cout << "Maximum Threads Per Block: " << deviceProp.maxThreadsPerBlock << std::endl;
    std::cout << "Shared Memory Per Block: " << deviceProp.sharedMemPerBlock << " bytes" << std::endl;

    return 0;
}
