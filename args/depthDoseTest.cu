#include <cuda_runtime.h>
#include <helper_cuda.h>

__global__ void
d_depthDose(float* output, int nPoints, cudaTextureObject_t texObj)
{
    uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
    float u = (float) x / nPoints;
    float voxel = tex1D<float>(texObj, u);
    output[x] = voxel;
}

extern "C"
void depthDose(dim3 gridSize, dim3 blockSize, float* output, \
    uint nPoints, cudaTextureObject_t& texObj)
{
    d_depthDose<<<gridSize, blockSize>>>(output, nPoints, texObj);
}