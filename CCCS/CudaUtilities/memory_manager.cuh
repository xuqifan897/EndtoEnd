#ifndef __MEMORY_MANAGER_CUH__
#define __MEMORY_MANAGER_CUH__
#include <stdlib.h>
#include <cmath>
#include <unordered_map>

#include <cuda.h>
#include <cuda_runtime.h>
#include "helper_cuda.h"

/* implements manager/wrapper around cuda device memory allocations with test functions and user customizable
 * memory usage limits
 */

#define KB 1024.
#define MB 1048576.
#define GB 1073741824.
#define TB 1099511627776.

class MemoryManager {
    protected:
        int64_t peak_mem_usage = 0;
        int64_t mem_usage = 0;
        int64_t mem_limit;
        void log_alloc(void* addr, int64_t usage);
        void log_dealloc(void* addr);

        std::unordered_map<void*, int64_t> alloc_table;

    public:
        MemoryManager() {}
        MemoryManager(int64_t mem_limit) : mem_limit{mem_limit} {}

        void setMemoryLimit(int64_t mem_limit) { this->mem_limit = mem_limit; };

        int64_t peakMemoryUsed() const {
            return mem_usage;
        }
        int64_t memoryUsed() const {
            return mem_usage;
        }
        int64_t memoryRemaining() const {
            return mem_limit - memoryUsed();
        }
        // test if available memory is greater than query
        int memoryAvailable(int64_t mem) const {
            return bool(memoryRemaining()>=mem);
        }


        //////////////////
        // ALLOC MEMORY //
        //////////////////
        template <typename T>
        void cudaMalloc(T **devPtr, size_t size) {
            cudaError_t err = ::cudaMalloc(devPtr, size);
            checkCudaErrors(err);
            log_alloc((void*)*devPtr, size);
        }

        void cudaMallocArray(cudaArray_t *array, const struct cudaChannelFormatDesc *desc, size_t width, size_t height=0, unsigned int flags=0);
        void cudaMalloc3DArray(struct cudaArray **array, const struct cudaChannelFormatDesc *desc, struct cudaExtent extent, unsigned int flags=0);
        void cudaMallocPitch(void **devPtr, size_t *pitch, size_t width, size_t height);
        void cudaMalloc3D(struct cudaPitchedPtr *pitchedDevPtr, struct cudaExtent extent);


        /////////////////
        // FREE MEMORY //
        /////////////////
        void cudaFree(void *devPtr);
        void cudaFreeArray(cudaArray_t array);
};

#endif // __MEMORY_MANAGER_CUH__
