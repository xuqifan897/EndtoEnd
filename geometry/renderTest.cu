#ifndef _RENDERTEST_KERNEL_CU_
#define _RENDERTEST_KERNEL_CU_

#include <helper_cuda.h>
#include <helper_math.h>
#include <assert.h>

__global__ void
d_render(float* output, int y_points, int z_points, float y_scale, \
    float z_scale, float x, cudaTextureObject_t texObj)
{
    uint y = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
    uint z = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;

    float u = y / (float) y_points * y_scale;
    float v = z / (float) z_points * z_scale;

    float voxel = tex3D<float>(texObj, v, u, x);
    uint i = __umul24(z, y_points) + y;
    output[i] = voxel;
}


extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, \
    float* output, int y_points, int z_points, float y_scale, float z_scale, \
    float x, cudaTextureObject_t& texObj)
{
    d_render<<<gridSize, blockSize>>>(output, y_points, z_points, \
        y_scale, z_scale, x, texObj);
}


__global__ void
d_convolve_kernel(float* convolvedFluenceMap, float* extendedFluenceMap, \
    float* convolutionKernel, int ConvRad, uint globalPitch, uint target_prepend, \
    uint source_prepend, uint kernel_prepend)
{   
    extern __shared__ float shared[];
    float* sharedResult = shared;
    float* sharedConvolutionKernel = &(shared[blockDim.x * blockDim.y]);
    float* extendSub = &(sharedConvolutionKernel[4*ConvRad*ConvRad]);
    assert(blockDim.x==ConvRad && blockDim.y==ConvRad);
    uint sharedPitch = 3 * ConvRad;
    // then extendSub is 9 times the block size
    uint global_row_base = blockIdx.x * ConvRad + source_prepend;
    uint global_col_base = blockIdx.y * ConvRad + source_prepend;
    
    for (int i=0; i<3; i++)
    {
        uint row = i * ConvRad + threadIdx.x;
        uint global_row = global_row_base + row;
        for (int j=0; j<3; j++)
        {
            uint col = j * ConvRad + threadIdx.y;
            uint global_col = global_col_base + col;
            size_t shared_idx = row * sharedPitch + col;
            size_t global_idx = global_row * globalPitch + global_col;
            extendSub[shared_idx] = extendedFluenceMap[global_idx];
        }
    }

    uint convolutionKernelPitch = 2*ConvRad;
    for (int i=0; i<2; i++)
    {
        uint row = i * ConvRad + threadIdx.x;
        for (int j=0; j<2; j++)
        {
            uint col = j * ConvRad + threadIdx.y;
            size_t idx = row * convolutionKernelPitch + col;
            sharedConvolutionKernel[idx] = convolutionKernel[idx];
        }
    }
    __syncthreads();

    uint sharedResultIdx = threadIdx.x * ConvRad + threadIdx.y;
    sharedResult[sharedResultIdx] = 0;
    for (int i=0; i<2*ConvRad-1; i++)
    {
        uint sharedRowIdx = threadIdx.x + i + 1;
        for (int j=0; j<2*ConvRad-1; j++)
        {
            uint sharedColIdx = threadIdx.y + j + 1;
            uint sharedIdx = sharedRowIdx * sharedPitch + sharedColIdx;
            uint sharedKernelIdx = (i + kernel_prepend) * convolutionKernelPitch + j + kernel_prepend;
            sharedResult[sharedResultIdx] += sharedConvolutionKernel[sharedKernelIdx] * extendSub[sharedIdx];
        }
    }
    // write back to global memory
    uint globalResultIdx = (blockIdx.x * blockDim.x + threadIdx.x + target_prepend) * (gridDim.y * \
        blockDim.y) + blockIdx.y * blockDim.y + threadIdx.y + target_prepend;
    convolvedFluenceMap[globalResultIdx] = sharedResult[sharedResultIdx];
}


__global__ void
d_convolve_kernel_no_shared(float* convolvedFluenceMap, float* extendedFluenceMap, \
    float* convolutionKernel, int ConvRad, uint globalPitch, uint target_prepend, \
    uint source_prepend, uint kernel_prepend)
{
    uint x = blockDim.x * blockIdx.x + threadIdx.x;
    uint y = blockDim.y * blockIdx.y + threadIdx.y;
    float value = 0;
    uint convolutionKernelPitch = 2 * ConvRad;
    for (int i=0; i<2*ConvRad-1; i++)
    {
        uint extended_x = x + 1 + i + source_prepend;
        for (int j=0; j<2*ConvRad-1; j++)
        {
            uint extended_y = y + 1 + j + source_prepend;
            uint extended_idx = extended_x * globalPitch + extended_y;
            uint convolutionKernel_idx = (i + kernel_prepend) * convolutionKernelPitch + j + kernel_prepend;
            value += extendedFluenceMap[extended_idx] * convolutionKernel[convolutionKernel_idx];
        }
    }
    uint result_pitch = blockDim.y * gridDim.y;
    uint result_idx = (x + target_prepend) * result_pitch + y + target_prepend;
    convolvedFluenceMap[result_idx] = value;
}


extern "C"
void convolve_kernel(dim3 gridSize, dim3 blockSize, cudaStream_t stream, float* convolvedFluenceMap, \
    float* extendedFluenceMap, float* convolutionKernel, int ConvRad, int globalPitch,  
    int target_prepend, int source_prepend, int kernel_prepend)
{
    uint sharedMemorySize = 14 * ConvRad * ConvRad * sizeof(float);
    d_convolve_kernel<<<gridSize, blockSize, sharedMemorySize, stream>>>(convolvedFluenceMap, \
        extendedFluenceMap, convolutionKernel, ConvRad, globalPitch, target_prepend, \
        source_prepend, kernel_prepend);

    // d_convolve_kernel_no_shared<<<gridSize, blockSize, 0, stream>>>(convolvedFluenceMap, \
    //     extendedFluenceMap, convolutionKernel, ConvRad, globalPitch, \
    //     target_prepend, source_prepend, kernel_prepend);
}

#endif