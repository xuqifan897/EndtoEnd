#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <assert.h>

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

__global__ void
d_convolution_kernel_init(float* data, int radius, int pitch, \
    float A, float B, float a, float b, float pixelSize, int n_samples_per_dim)
{
    int x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
    float centerCoord_x = x + 1 - radius;
    float centerCoord_y = y + 1 - radius;

    float reference_x = (float)(n_samples_per_dim - 1) / 2;
    float reference_y = reference_x;
    float avg = 0;
    for (int i=0; i<n_samples_per_dim; i++)
    {
        float coord_x = centerCoord_x + (float)(i - reference_x)/n_samples_per_dim;
        for (int j=0; j<n_samples_per_dim; j++)
        {
            float coord_y = centerCoord_y + (float)(j - reference_y)/n_samples_per_dim;
            float rho = sqrt(coord_x * coord_x + coord_y * coord_y) * pixelSize;
            avg += (exp(-a*rho)*A + exp(-b*rho)*B)/rho;
        }
    }
    avg /= (n_samples_per_dim * n_samples_per_dim);

    uint idx = x * pitch + y;
    data[idx] = avg;
}

__global__ void
d_convolution_kernel_normalization(float* data, int coord0, int coord1, int pitch)
{
    int x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
    int reference_idx = coord0 * pitch + coord1;
    float reference_value = data[reference_idx];
    int idx = x * pitch + y;
    float value = data[idx];
    value /= reference_value;
    data[idx] = value;
}

__global__ void
d_convolution_kernel_init_even(float* data, int radius, int pitch, \
    float A, float B, float a, float b, float pixelSize)
{
    float center = ((float)pitch - 1) / 2;
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    uint idx = idx_x * pitch + idx_y;
    
    float centerCoord_x = ((float)idx_x - center) * pixelSize;
    float centerCoord_y = ((float)idx_y - center) * pixelSize;
    float rho = sqrt(centerCoord_x * centerCoord_x + centerCoord_y * centerCoord_y);
    float value = A * exp(-a * rho) + B * exp(-b * rho);
    data[idx] = value;
}

extern "C"
void convolution_kernel_init(dim3 gridSize, dim3 blockSize, float* data, int radius, \
    int pitch, float A, float B, float a, float b, float pixelSize, int n_samples_per_dim)
{   
    // assert(n_samples_per_dim % 2 == 0);
    // d_convolution_kernel_init<<<gridSize, blockSize>>>( \
    //     data, radius, pitch, A, B, a, b, pixelSize, n_samples_per_dim);

    d_convolution_kernel_init_even<<<gridSize, blockSize>>>( \
        data, radius, pitch, A, B, a, b, pixelSize);
    d_convolution_kernel_normalization<<<gridSize, blockSize>>>( \
        data, radius-1, radius-1, pitch);
}