#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"

using namespace E2E;

__global__ void
d_doseSum(float* result, float** sources, uint num_beams, uint size)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx > size)
        return;
    
    float value = 0;
    for (uint i=0; i<num_beams; i++)
        value += sources[i][idx];
    result[idx] = value;
}

extern "C"
void doseSum(float* result, float** sources, uint num_beams, uint size, cudaStream_t stream)
{
    // size is the size of each source
    dim3 blockSize(64);
    dim3 gridSize(ceil((float)size / blockSize.x));
    d_doseSum<<<gridSize, blockSize, 0, stream>>>(result, sources, num_beams, size);
}


__device__ void warpReduce(volatile float* sdata, uint tid)
{
    if (REDUCTION_BLOCK_SIZE >= 64)
        sdata[tid] += sdata[tid + 32];
    if (REDUCTION_BLOCK_SIZE >= 32)
        sdata[tid] += sdata[tid + 16];
    if (REDUCTION_BLOCK_SIZE >= 16)
        sdata[tid] += sdata[tid + 8];
    if (REDUCTION_BLOCK_SIZE >= 8)
        sdata[tid] += sdata[tid + 4];
    if (REDUCTION_BLOCK_SIZE >= 4)
        sdata[tid] += sdata[tid + 2];
    if (REDUCTION_BLOCK_SIZE >= 2)
        sdata[tid] += sdata[tid + 1];
}


__global__ void
d_Reduction(float* out, float* in, uint size)
{
    __shared__ float sdata[REDUCTION_BLOCK_SIZE];
    uint tid = threadIdx.x;
    uint i = blockIdx.x * REDUCTION_BLOCK_SIZE + tid;
    uint gridSize = REDUCTION_BLOCK_SIZE * gridDim.x;
    sdata[tid] = 0;

    while (i < size)
    {
        sdata[tid] += in[i];
        i += gridSize;
    }
    __syncthreads();

    if (REDUCTION_BLOCK_SIZE >= 512)
    {
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
    }

    if (REDUCTION_BLOCK_SIZE >= 256)
    {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
    }

    if (REDUCTION_BLOCK_SIZE >= 128)
    {
        if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
    }

    if (tid < 32)
        warpReduce(sdata, tid);
    if (tid == 0)
        out[blockIdx.x] = sdata[0];
}

extern "C"
void Reduction(dim3 gridSize, dim3 blockSize, float* out, \
    float* source, uint size, cudaStream_t stream)
{
    d_Reduction<<<gridSize, blockSize, 0, stream>>>(out, source, size);
}