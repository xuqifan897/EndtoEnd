#ifndef _RENDERTEST_KERNEL_CU_
#define _RENDERTEST_KERNEL_CU_

#include <helper_cuda.h>
#include <helper_math.h>

__global__ void
d_render(float* output, int y_points, int z_points, float x, cudaTextureObject_t texObj)
{
    uint y = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
    uint z = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;

    float u = y / (float) y_points;
    float v = z / (float) z_points;

    float voxel = tex3D<float>(texObj, v, u, x);
    {
        uint i = __umul24(z, y_points) + y;
        output[i] = voxel;
    }
}


extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, \
    float* output, int y_points, int z_points, float x, cudaTextureObject_t& texObj)
{
    d_render<<<gridSize, blockSize>>>(output, y_points, z_points, x, texObj);
}
#endif