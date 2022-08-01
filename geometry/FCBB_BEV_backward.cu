#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"

using namespace E2E;


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
d_BEVDoseBackward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    float phantom_size0, float phantom_size1, float phantom_size2, \
    float phantom_iso0, float phantom_iso1, float phantom_iso2, \
    uint fluence_map_dimension, \
    float* d_convolved_fluence_grad, \
    cudaTextureObject_t phantom_texture, \
    float* d_FCBB_BEV_dose_grad, \
    float max_depth, cudaTextureObject_t depthDoseTexture)
{
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (idx_x >= fluence_map_dimension || idx_y >= fluence_map_dimension) return;

    float sampling_step = (sampling_range_end - sampling_range_start) / (sampling_points - 1);
    float x_from_center = ((float)idx_x - (float)(fluence_map_dimension - 2) / 2) * pixel_size;
    float y_from_center = ((float)idx_y - (float)(fluence_map_dimension - 2) / 2) * pixel_size;
    float r0 = sqrt(SAD*SAD + x_from_center*x_from_center +  y_from_center*y_from_center);
    float radiological_path_step = r0 * sampling_step / SAD;
    float radiological_path = 0;

    float BEV_source_to_isocenter[3]{0, 0, -SAD};
    float PVCS_coords[3];
    d_BEV_to_PVCS(PVCS_coords, BEV_source_to_isocenter, zenith, azimuth);
    float PVCS_source[3]{phantom_iso0 + PVCS_coords[0], phantom_iso1 + PVCS_coords[1], \
        phantom_iso2 + PVCS_coords[2]};
    float BEV_reference[3]{x_from_center, y_from_center, SAD};
    float PVCS_reference[3];
    d_BEV_to_PVCS(PVCS_reference, BEV_reference, zenith, azimuth);

    uint FCBB_BEV_dose_grad_pitch = fluence_map_dimension * fluence_map_dimension;
    uint FCBB_BEV_dose_grad_offset = idx_x * fluence_map_dimension + idx_y;
    float fluence_grad = 0;
    for (uint i=0; i<sampling_points; i++)
    {
        float length = sampling_range_start + i * sampling_step;
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
        float BEV_dose_grad = d_FCBB_BEV_dose_grad[i * FCBB_BEV_dose_grad_pitch + FCBB_BEV_dose_grad_offset];
        fluence_grad += tex1D<float>(depthDoseTexture, normalized_radiological_path) / (scale * scale) * BEV_dose_grad;
    }
    d_convolved_fluence_grad[FCBB_BEV_dose_grad_offset] = fluence_grad * \
        !(idx_x == fluence_map_dimension-1 || idx_y == fluence_map_dimension-1);
    // the last pixel of the convolved fluence map is discarded
}

extern "C"
void BEVDoseBackward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    float phantom_size[3], float phantom_iso[3], \
    float* d_convolved_fluence_grad, \
    cudaTextureObject_t phantom_texture, \
    float* d_FCBB_BEV_dose_grad, \
    FCBBkernel* FCBB_kernel, \
    cudaStream_t stream)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    dim3 blockSize(16, 16);
    dim3 gridSize(convolved_fluence_map_dimension / blockSize.x, \
        convolved_fluence_map_dimension / blockSize.y);
    d_BEVDoseBackward<<<gridSize, blockSize, 0, stream>>>(zenith, azimuth, SAD, pixel_size, \
        sampling_range_start, sampling_range_end, sampling_points, \
        phantom_size[0], phantom_size[1], phantom_size[2], \
        phantom_iso[0], phantom_iso[1], phantom_iso[2], \
        convolved_fluence_map_dimension, \
        d_convolved_fluence_grad, \
        phantom_texture, \
        d_FCBB_BEV_dose_grad, \
        (*FCBB_kernel).max_depth, (*FCBB_kernel).tex);
}