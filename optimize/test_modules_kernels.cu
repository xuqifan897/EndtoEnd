#include <cuda_runtime.h>
#include <helper_cuda.h>


__global__ void
d_moduleTestHostTrilinear(cudaTextureObject_t tex, float* d_result, uint dim0, uint dim1, \
    float level2, float shift0, float shift1)
{
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (idx_x >= dim0 || idx_y >= dim1)
        return;
    uint idx = idx_x * dim1 + idx_y;

    float x_norm = (float)idx_x / dim0 + shift0;
    float y_norm = (float)idx_y / dim1 + shift1;
    d_result[idx] = tex3D<float>(tex, level2, y_norm, x_norm);
}


extern "C"
void moduleTestHostTrilinear(cudaTextureObject_t tex, float* d_result, uint dim0, uint dim1, \
    float level2, float shift0, float shift1)
{
    dim3 blockSize(16, 16);
    dim3 gridSize(ceil((float)dim0 / blockSize.x), ceil((float)dim1 / blockSize.y));
    d_moduleTestHostTrilinear<<<blockSize, gridSize>>>(tex, d_result, dim0, dim1, \
        level2, shift0, shift1);
}


__global__ void
d_moduleTestHostLinear(cudaTextureObject_t tex, float* d_result, uint dim, float shift)
{
    uint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= dim) return;
    float idx_normalized = (float)idx / dim + shift;
    d_result[idx] = tex1D<float>(tex, idx_normalized);
}


extern "C"
void moduleTestHostLinear(cudaTextureObject_t tex, float* d_result, uint dim, float shift)
{
    dim3 blockSize(32);
    dim3 gridSize(ceil((float)dim / blockSize.x));
    d_moduleTestHostLinear<<<gridSize, blockSize>>>(tex, d_result, dim, shift);
}