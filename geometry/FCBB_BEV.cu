#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <assert.h>
#include "geom.h"
#include "args.h"

using namespace E2E;
using namespace std;

// for debug purposes
float* E2E::HU_debug = nullptr;
float* E2E::dose_debug = nullptr;


__device__ void inline
d_BEV_to_PVCS(float* PVCS_coords, float* BEV_coords, float zenith, float azimuth)
{
    float sin_zenith = sin(zenith);
    float cos_zenith = cos(zenith);
    float sin_azimuth = sin(azimuth);
    float cos_azimuth = cos(azimuth);

    // // left hand coordinates
    // float x_hat[3]{cos_zenith * cos_azimuth, - cos_zenith * sin_azimuth, sin_zenith};
    // float y_hat[3]{sin_azimuth, cos_azimuth, 0};
    // float z_hat[3]{- sin_zenith * cos_azimuth, sin_zenith * sin_azimuth, cos_zenith};

    // right hand coordinates
    float x_hat[3]{cos_zenith * cos_azimuth, cos_zenith * sin_azimuth, -sin_zenith};
    float y_hat[3]{-sin_azimuth, cos_azimuth, 0};
    float z_hat[3]{sin_zenith * cos_azimuth, sin_zenith * sin_azimuth, cos_zenith};

    // store temporary value
    PVCS_coords[0] = BEV_coords[0] * x_hat[0] + BEV_coords[1] * y_hat[0] + BEV_coords[2] * z_hat[0];
    PVCS_coords[1] = BEV_coords[0] * x_hat[1] + BEV_coords[1] * y_hat[1] + BEV_coords[2] * z_hat[1];
    PVCS_coords[2] = BEV_coords[0] * x_hat[2] + BEV_coords[1] * y_hat[2] + BEV_coords[2] * z_hat[2];
}


__global__ void
d_BEVDoseForward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_step, uint sampling_points, \
    cudaSurfaceObject_t dose_surface, \
    float phantom_size0, float phantom_size1, float phantom_size2, \
    float phantom_iso0, float phantom_iso1, float phantom_iso2, \
    cudaTextureObject_t phantom_texture, \
    float max_depth, cudaTextureObject_t depthDose_texture, \
    float* convolved_fluence_map)
    // // for debug purposes
    // float* d_HU_debug, float* d_dose_debug)
{
    uint x_idx = blockDim.x * blockIdx.x + threadIdx.x;
    uint y_idx = blockDim.y * blockIdx.y + threadIdx.y;
    uint fluence_map_dimension = blockDim.x * gridDim.x;
    assert(fluence_map_dimension == blockDim.y * gridDim.y);

    float fluence = convolved_fluence_map[x_idx * fluence_map_dimension + y_idx];
    float x_from_center = ((float)x_idx - ((float)fluence_map_dimension - 1)/2) * pixel_size;
    float y_from_center = ((float)y_idx - ((float)fluence_map_dimension - 1)/2) * pixel_size;
    float radiological_path_step = sqrt(SAD*SAD + x_from_center*x_from_center + \
        y_from_center*y_from_center) * sampling_step / SAD;
    float radiological_path = 0;

    float BEV_source_to_isocenter[3]{0, 0, -SAD};
    float PVCS_coords[3];
    d_BEV_to_PVCS(PVCS_coords, BEV_source_to_isocenter, zenith, azimuth);
    float PVCS_source[3]{phantom_iso0 + PVCS_coords[0], phantom_iso1 + PVCS_coords[1], \
        phantom_iso2 + PVCS_coords[2]};
    float BEV_reference[3]{x_from_center, y_from_center, SAD};
    float PVCS_reference[3];
    d_BEV_to_PVCS(PVCS_reference, BEV_reference, zenith, azimuth);

    
    for (uint i=0; i<sampling_points; i++)
    {
        float length = i * sampling_step + sampling_range_start;
        float scale = length / SAD;
        PVCS_coords[0] = PVCS_source[0] + PVCS_reference[0] * scale;
        PVCS_coords[1] = PVCS_source[1] + PVCS_reference[1] * scale;
        PVCS_coords[2] = PVCS_source[2] + PVCS_reference[2] * scale;

        PVCS_coords[0] = PVCS_coords[0] / phantom_size0;
        PVCS_coords[1] = PVCS_coords[1] / phantom_size1;
        PVCS_coords[2] = PVCS_coords[2] / phantom_size2;

        // patient texture index follow (z, y, x) order
        float HU = tex3D<float>(phantom_texture, PVCS_coords[2], PVCS_coords[1], PVCS_coords[0]);
        radiological_path += HU * radiological_path_step;
        float normalized_radiological_path = radiological_path / max_depth;
        float dose = tex1D<float>(depthDose_texture, normalized_radiological_path) * \
            fluence * HU / (scale * scale);

        // // for debug purposes
        // uint debug_idx = (x_idx * fluence_map_dimension + y_idx) * sampling_points + i;
        // d_HU_debug[debug_idx] = HU;
        // d_dose_debug[debug_idx] = dose;

        // BEV_dose_surface follow (y, x, z) order
        surf3Dwrite(dose, dose_surface, y_idx*sizeof(float), x_idx, i, cudaBoundaryModeZero);
    }
}

extern "C"
void BEVDoseForward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    cudaSurfaceObject_t dose_surface, \
    float phantom_size[3], float phantom_iso[3], \
    cudaTextureObject_t phantom_texture, \
    float* convolved_fluence_map, \
    cudaStream_t stream=0)
{
    uint blockS = 16;
    dim3 blockSize(blockS, blockS);
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    float sampling_step = (sampling_range_end - sampling_range_start) / (sampling_points - 1);
    dim3 gridSize(convolved_fluence_map_dimension / blockS, \
        convolved_fluence_map_dimension / blockS);
    d_BEVDoseForward<<<gridSize, blockSize, 0, stream>>>( \
        zenith, azimuth, SAD, pixel_size, sampling_range_start, sampling_step, \
        sampling_points, dose_surface, \
        phantom_size[0], phantom_size[1], phantom_size[2], \
        phantom_iso[0], phantom_iso[1], phantom_iso[2], \
        phantom_texture, \
        (*FCBB6MeV).max_depth, (*FCBB6MeV).tex, \
        convolved_fluence_map);
        // // for debug purposes
        // HU_debug, dose_debug \
        );
}


__global__ void
d_writeSurface(cudaSurfaceObject_t surface, float* data)
{
    uint x = blockDim.x * blockIdx.x + threadIdx.x;
    uint y = blockDim.y * blockIdx.y + threadIdx.y;
    uint z = blockDim.z * blockIdx.z + threadIdx.z;

    uint y_pitch = blockDim.z * gridDim.z;
    uint x_pitch = blockDim.y * gridDim.y;
    uint data_idx = (x * x_pitch + y) * y_pitch + z;
    surf3Dwrite(data[data_idx], surface, z * sizeof(float), y, x, cudaBoundaryModeZero);
}


extern "C"
void writeSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data)
{
    d_writeSurface<<<gridSize, blockSize>>>(surface, data);
}


__global__ void
d_readTexture(cudaTextureObject_t texture, float* data)
{
    uint x = blockDim.x * blockIdx.x + threadIdx.x;
    uint y = blockDim.y * blockIdx.y + threadIdx.y;
    uint z = blockDim.z * blockIdx.z + threadIdx.z;

    uint x_dim = blockDim.x * gridDim.x;
    uint y_dim = blockDim.y * gridDim.y;
    uint z_dim = blockDim.z * gridDim.z;

    uint data_idx = (x * y_dim + y) * z_dim + z;

    float x_norm = (float)x / x_dim;
    float y_norm = (float)y / y_dim;
    float z_norm = (float)z / z_dim;

    data[data_idx] = tex3D<float>(texture, z_norm, y_norm, x_norm);
}

extern "C"
void readTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* data)
{
    d_readTexture<<<gridSize, blockSize>>>(texture, data);
}


__global__ void
d_readSurface(cudaSurfaceObject_t surface, float* data)
{
    /* for BEV dose surface, the cudaExtent order: (y, x, z) 
    logical order: (z, x, y) */
    uint x = blockDim.x * blockIdx.x + threadIdx.x;
    uint y = blockDim.y * blockIdx.y + threadIdx.y;
    uint z = blockDim.z * blockIdx.z + threadIdx.z;

    uint x_dim = blockDim.x * gridDim.x;
    uint y_dim = blockDim.y * gridDim.y;
    // uint z_dim = blockDim.z * gridDim.z;

    uint data_idx = (z * x_dim + x) * y_dim + y;
    // data[data_idx] = surf3Dread(surface, y*sizeof(float), x, z, cudaBoundaryModeTrap);
    surf3Dread(&(data[data_idx]), surface, y*sizeof(float), x, z);
}

extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data)
{
    d_readSurface<<<gridSize, blockSize>>>(surface, data);
}