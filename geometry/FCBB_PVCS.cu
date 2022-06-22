#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

using namespace std;
const float half = 0.5;

__device__ void inline
d_PVCS_to_BEV(float* BEV_coords, float* PVCS_coords, float zenith, float azimuth)
{
    float sin_zenith = sin(zenith);
    float cos_zenith = cos(zenith);
    float sin_azimuth = sin(azimuth);
    float cos_azimuth = cos(azimuth);

    float x_hat[3]{cos_zenith * cos_azimuth, -sin_azimuth, sin_zenith * cos_azimuth};
    float y_hat[3]{cos_zenith * sin_azimuth, cos_azimuth, sin_zenith * sin_azimuth};
    float z_hat[3]{-sin_zenith, 0, cos_zenith};

    BEV_coords[0] = PVCS_coords[0] * x_hat[0] + PVCS_coords[1] * y_hat[0] + PVCS_coords[2] * z_hat[0];
    BEV_coords[1] = PVCS_coords[0] * x_hat[1] + PVCS_coords[1] * y_hat[1] + PVCS_coords[2] * z_hat[1];
    BEV_coords[2] = PVCS_coords[0] * x_hat[2] + PVCS_coords[1] * y_hat[2] + PVCS_coords[2] * z_hat[2];
}


__global__ void
d_PVCSDoseForward(float voxel_size, \
    uint phantom_dim0, uint phantom_dim1, uint phantom_dim2, uint phantom_pitch, \
    float phantom_iso0, float phantom_iso1, float phantom_iso2, \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, float sampling_step, \
    float fluence_map_size_x, float fluence_map_size_y, \
    float* FCBB_PVCS_dose, \
    cudaTextureObject_t BEV_dose_texture)
{
    uint x_idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint y_idx = blockIdx.y * blockDim.y + threadIdx.y;
    uint z_idx = blockIdx.z * blockDim.z + threadIdx.z;

    if (x_idx >= phantom_dim0 || y_idx >= phantom_dim1 || z_idx >= phantom_dim2)
        return;

    // here, we assume the points locate at the center of each voxel in the phantom
    float PVCS_point[3]{((float)x_idx + half) * voxel_size, \
        ((float)y_idx + half) * voxel_size, ((float)z_idx + half) * voxel_size};
    float PVCS_point_minus_iso[3]{PVCS_point[0] - phantom_iso0, \
        PVCS_point[1] - phantom_iso1, PVCS_point[2] - phantom_iso2};

    float BEV_point_minus_iso[3];
    d_PVCS_to_BEV(BEV_point_minus_iso, PVCS_point_minus_iso, zenith, azimuth);

    float BEV_point[3]{BEV_point_minus_iso[0], BEV_point_minus_iso[1], BEV_point_minus_iso[2] + SAD};

    float fluence_map_size_x_at_this_depth = fluence_map_size_x * BEV_point[2] / SAD;
    float fluence_map_size_y_at_this_depth = fluence_map_size_y * BEV_point[2] / SAD;
    // float normalized_BEV_point[3]{BEV_point[0] / fluence_map_size_x_at_this_depth + half, \
    //     BEV_point[1] / fluence_map_size_y_at_this_depth + half, \
    //     (BEV_point[2] + half * sampling_step) / (sampling_end - sampling_start + sampling_step)};

    float normalized_BEV_point[3]{BEV_point[0] / fluence_map_size_x_at_this_depth + half, \
        BEV_point[1] / fluence_map_size_y_at_this_depth + half, \
        (BEV_point[2] - sampling_start + half * sampling_step) / \
        (sampling_end - sampling_start + sampling_step)};
    
    // BEV_dose logical order: (z, x, y), cudaExtent order: (y, x, z)
    float value = tex3D<float>(BEV_dose_texture, normalized_BEV_point[1], \
        normalized_BEV_point[0], normalized_BEV_point[2]);

    uint idx = (x_idx * phantom_dim1 + y_idx) * phantom_pitch + z_idx;
    FCBB_PVCS_dose[idx] = value;
}


extern "C" void
PVCSDoseForward(float voxel_size, uint phantom_dim[3], uint phantom_pitch, \
    float phantom_iso[3], \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_size_x, float fluence_map_size_y, 
    float* FCBB_PVCS_dose, \
    cudaTextureObject_t BEV_dose_texture, \
    cudaStream_t stream=0)
{
    // the logical order of BEV_dose_texture is (z, x, y), cudaExtent order (y, x, z)
    float sampling_step = (sampling_end - sampling_start) / (sampling_points - 1);
    dim3 blockSize(1, 16, 16);
    dim3 gridSize;
    gridSize.x = ceil((float)(phantom_dim[0]) / blockSize.x);
    gridSize.y = ceil((float)(phantom_dim[1]) / blockSize.y);
    gridSize.z = ceil((float)(phantom_dim[2]) / blockSize.z);

    d_PVCSDoseForward<<<gridSize, blockSize, 0, stream>>>(voxel_size, \
        phantom_dim[0], phantom_dim[1], phantom_dim[2], phantom_pitch, \
        phantom_iso[0], phantom_iso[1], phantom_iso[2], \
        zenith, azimuth, SAD, \
        sampling_start, sampling_end, sampling_step, \
        fluence_map_size_x, fluence_map_size_y, \
        FCBB_PVCS_dose, BEV_dose_texture);
}


__global__ void
d_testWritePVCSSurface(cudaSurfaceObject_t surface)
{
    uint x_idx = blockDim.x * blockIdx.x + threadIdx.x;
    uint y_idx = blockDim.y * blockIdx.y + threadIdx.y;
    uint z_idx = blockDim.z * blockIdx.z + threadIdx.z;

    float value = (float)x_idx + y_idx - z_idx;
    // cudaExtent order: (z, y, x)
    surf3Dwrite(value, surface, z_idx * sizeof(float), y_idx, x_idx, cudaBoundaryModeZero);
}

extern "C"
void testWritePVCSSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface)
{
    d_testWritePVCSSurface<<<gridSize, blockSize>>>(surface);
}

__global__ void
d_testReadPVCSTexture(cudaTextureObject_t texture, float* output)
{
    uint x_idx = blockDim.x * blockIdx.x + threadIdx.x;
    uint y_idx = blockDim.y * blockIdx.y + threadIdx.y;
    uint z_idx = blockDim.z * blockIdx.z + threadIdx.z;

    // uint x_dim = blockDim.x * gridDim.x;
    uint y_dim = blockDim.y * gridDim.y;
    uint z_dim = blockDim.z * gridDim.z;

    // cudaExtent order: (z, y, x)
    float value = tex3D<float>(texture, z_idx, y_idx, x_idx);
    uint output_idx = (x_idx * y_dim + y_idx) * z_dim + z_idx;
    output[output_idx] = value;
}

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output)
{
    d_testReadPVCSTexture<<<gridSize, blockSize>>>(texture, output);
}