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
d_Reduction(float* out, float* in, uint size, uint idx)
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
        out[blockIdx.x + idx] = sdata[0];
}

extern "C"
void Reduction(dim3 gridSize, dim3 blockSize, float* out, \
    float* source, uint size, uint idx, cudaStream_t stream)
{
    d_Reduction<<<gridSize, blockSize, 0, stream>>>(out, source, size, idx);
}

__global__ void
d_elementWiseSquare(float* d_output, float* d_input, uint size)
{
    uint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx>=size) return;
    d_output[idx] = d_input[idx] * d_input[idx];
}

extern "C"
void elementWiseSquare(dim3 gridSize, dim3 blockSize, float* d_output, \
    float* d_input, uint size, cudaStream_t stream)
{
    d_elementWiseSquare<<<gridSize, blockSize, 0, stream>>>(d_output, d_input, size);
}

__global__ void
d_fluenceMapUpdate(uint fluence_map_dimension, uint extended_fluence_map_dimension, uint extended_fluence_map_offset, \
    uint idx, float* d_extended_fluence_map, float* d_fluence_map_grad, float* d_norm, float step_size)
{
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (idx_x >= fluence_map_dimension || idx_y >= fluence_map_dimension)
        return;
    
    float real_norm = d_norm[idx];
    real_norm = real_norm / (fluence_map_dimension * fluence_map_dimension);
    real_norm = sqrt(real_norm);

    uint fluence_map_idx = idx_x * fluence_map_dimension + idx_y;
    uint extended_fluence_map_idx = (idx_x + extended_fluence_map_offset) * \
        extended_fluence_map_dimension + idx_y + extended_fluence_map_offset;
    d_extended_fluence_map[extended_fluence_map_idx] -= d_fluence_map_grad[fluence_map_idx] * step_size / real_norm;
}

extern "C"
void fluenceMapUpdate(uint fluence_map_dimension, uint extended_fluence_map_dimension, uint extended_fluence_map_offset, \
    uint idx, float* d_extended_fluence_map, float* d_fluence_map_grad, float* d_norm, float step_size, \
    cudaStream_t stream)
{
    dim3 blockSize(32, 32);
    dim3 gridSize((uint)ceil((float)fluence_map_dimension / blockSize.x), \
        (uint)ceil((float)fluence_map_dimension / blockSize.y));
    d_fluenceMapUpdate<<<gridSize, blockSize, 0, stream>>>(fluence_map_dimension, extended_fluence_map_dimension, extended_fluence_map_offset, \
        idx, d_extended_fluence_map, d_fluence_map_grad, d_norm, step_size);
}