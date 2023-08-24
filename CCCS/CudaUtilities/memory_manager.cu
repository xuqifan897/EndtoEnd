#include "memory_manager.cuh"

#include <cstdio>
#include <exception.h>
#include <cstring>

void MemoryManager::log_alloc(void* addr, int64_t usage) {
    if (!memoryAvailable(usage)) {
        char msg[300];
        sprintf(msg, "Request to allocate %d bytes of device memory exceeded remaining memory (%d bytes)", usage, memoryRemaining());
        throw std::runtime_error(msg);
    }
    mem_usage += usage;
    peak_mem_usage = max(peak_mem_usage, mem_usage);
    alloc_table.emplace(addr, usage);
    /* printf("+++ added memory usage to table: %p (%d)\n", addr, usage); */
}
void MemoryManager::log_dealloc(void* addr) {
    int64_t usage;
    /* printf("=== trying to locate memory usage in table: %p\n", addr); */
    try {
        usage = alloc_table.at(addr);
    } catch (const std::out_of_range& e) {
        printf("Request to free GPU memory failed because record was absent from allocation table.");
        throw e;
    }
    mem_usage -= usage;
    alloc_table.erase(addr);
    /* printf("--- removed memory usage from table: %p (%d)\n", addr, usage); */
}

void MemoryManager::cudaMallocArray(struct cudaArray **array, const struct cudaChannelFormatDesc *desc, size_t width, size_t height, unsigned int flags) {
    int elementsize = (desc->w+desc->x+desc->y+desc->z)/8;
    cudaError_t err = ::cudaMallocArray(array, desc, width, height, flags);
    checkCudaErrors(err);
    log_alloc((void*)*array, width*max(1,(int)height)*elementsize);
}

void MemoryManager::cudaMalloc3DArray(struct cudaArray **array, const struct cudaChannelFormatDesc *desc, struct cudaExtent extent, unsigned int flags) {
    int elementsize = (desc->w+desc->x+desc->y+desc->z)/8;
    cudaError_t err = ::cudaMalloc3DArray(array, desc, extent, flags);
    checkCudaErrors(err);
    log_alloc((void*)*array, extent.width*max(1,(int)extent.height)*max(1,(int)extent.depth)*elementsize);
}

void MemoryManager::cudaMallocPitch(void **devPtr, size_t *pitch, size_t width, size_t height) {
    cudaError_t err = ::cudaMallocPitch(devPtr, pitch, width, height);
    checkCudaErrors(err);
    log_alloc((void*)*devPtr, *pitch*height);
}

void MemoryManager::cudaMalloc3D(struct cudaPitchedPtr *pitchedDevPtr, struct cudaExtent extent) {
    cudaError_t err = ::cudaMalloc3D(pitchedDevPtr, extent);
    checkCudaErrors(err);
    log_alloc((void*)pitchedDevPtr, pitchedDevPtr->pitch*extent.height*extent.depth);
}


void MemoryManager::cudaFree(void *devPtr) {
    cudaError_t err = ::cudaFree(devPtr);
    checkCudaErrors(err);
    log_dealloc(devPtr);
}
void MemoryManager::cudaFreeArray(struct cudaArray *array) {
    cudaError_t err = ::cudaFreeArray(array);
    checkCudaErrors(err);
    log_dealloc((void*)array);
}
