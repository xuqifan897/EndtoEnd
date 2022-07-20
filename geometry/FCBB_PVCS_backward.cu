#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>
#include "geom.h"

using namespace E2E;
using namespace std;

bool* E2E::valid_debug = nullptr;

const float half = 0.5;

__global__ void
d_FCBBPVCSDoseGrad(cudaSurfaceObject_t FCBB_PVCS_dose_grad_surface, float* elementWiseLoss, \
    float* d_FCBB_PVCS_dose, float* PTV_weight, float* PTV_target, \
    float* OAR_weight, float* OAR_target, uint dimension0, uint dimension1, uint dimension2)
    // // for debug purposes
    // float* debug_result)
{
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    uint idx_z = blockIdx.z * blockDim.z + threadIdx.z;

    size_t idx = (idx_x * dimension1 + idx_y) * dimension2 + idx_z;
    if ((idx_x >= dimension0) || (idx_y >= dimension1) || (idx_z >= dimension2))
        return;
    float PTV_diff = d_FCBB_PVCS_dose[idx] - PTV_target[idx];
    float OAR_diff = d_FCBB_PVCS_dose[idx] - OAR_target[idx];
    float grad_value = 2 * PTV_weight[idx] * PTV_diff + 2 * OAR_weight[idx] * OAR_diff * (OAR_diff > 0);
    elementWiseLoss[idx] = PTV_weight[idx] * PTV_diff * PTV_diff + \
        OAR_weight[idx] * OAR_diff * OAR_diff * (OAR_diff > 0);

    // // for debug purposes
    // debug_result[idx] = grad_value;

    surf3Dwrite(grad_value, FCBB_PVCS_dose_grad_surface, idx_z*sizeof(float), \
        idx_y, idx_x, cudaBoundaryModeTrap);
}


extern "C"
void FCBBPVCSDoseGrad(cudaSurfaceObject_t FCBB_PVCS_dose_grad_surface, float* elementWiseLoss, \
    float* d_FCBB_PVCS_dose, float* PTV_weight, float* PTV_target, \
    float* OAR_weight, float* OAR_target, uint dimension[3], cudaStream_t stream)
{
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(ceil((float)dimension[0] / blockSize.x), \
        ceil((float)dimension[1] / blockSize.y), ceil((float)dimension[2] / blockSize.z));
    d_FCBBPVCSDoseGrad<<<gridSize, blockSize, 0, stream>>>( \
        FCBB_PVCS_dose_grad_surface, elementWiseLoss, \
        d_FCBB_PVCS_dose, PTV_weight, PTV_target, OAR_weight, OAR_target, \
        dimension[0], dimension[1], dimension[2]);
        // // for debug purposes
        // dose_debug);
}

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
d_PVCSDoseBackward(float voxel_size, float phantom_iso0, float phantom_iso1, float phantom_iso2, \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_pixel_size, uint fluence_map_dimension, \
    float* d_FCBB_BEV_dose_grad, cudaTextureObject_t PVCSDoseGradTexture)
{
    float sampling_step = (sampling_end - sampling_start) / (sampling_points - 1);
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    uint idx_z = blockIdx.z * blockDim.z + threadIdx.z;
    
    /* first, calculate the BEV coordinates of the BEV grid points. 
    The same as PVCS_forward, the the dimension of the fluence map here is 
    FM_dimension + 2 * FM_convolution_radius, with the last pixel discarded.
    So the effective index range is 0, 1, ..., FM_dimension + 2 * FM_convolution_radius - 2 */
    float projection_z = (sampling_start + sampling_step * idx_z);
    float fluence_map_pixel_size_at_this_depth = fluence_map_pixel_size * projection_z / SAD;
    float middle_point = (float)(fluence_map_dimension - 2) / 2;
    float BEV_coordinates_of_BEV_grid_points[3]{ \
        ((float)idx_x - middle_point) * fluence_map_pixel_size_at_this_depth, \
        ((float)idx_y - middle_point) * fluence_map_pixel_size_at_this_depth, \
        projection_z - SAD};
    
    // second, we convert the BEV grid points to PVCS coordinates
    float PVCS_coordinates_of_BEV_grid_points[3];
    d_BEV_to_PVCS(PVCS_coordinates_of_BEV_grid_points, BEV_coordinates_of_BEV_grid_points, zenith, azimuth);
    PVCS_coordinates_of_BEV_grid_points[0] += phantom_iso0;
    PVCS_coordinates_of_BEV_grid_points[1] += phantom_iso1;
    PVCS_coordinates_of_BEV_grid_points[2] += phantom_iso2;

    // third, we get the coordinates and grad values of the 8 neighboring PVCS grid points
    // at the beginning, we convert the physical coordinates of 
    // PVCS_coordinates_of_BEV_grid_points to unitless values
    int PVCS_index_of_BEV_grid_points[3]{(int)floorf(PVCS_coordinates_of_BEV_grid_points[0] / voxel_size - half), \
        (int)floorf(PVCS_coordinates_of_BEV_grid_points[1] / voxel_size - half), \
        (int)floorf(PVCS_coordinates_of_BEV_grid_points[2] / voxel_size - half)};
    int temp_PVCS_index_of_PVCS_grid_points[3];
    float PVCS_coordinates_of_PVCS_grid_points[3];
    float BEV_coordinates_of_PVCS_grid_points[3];
    float BEV_index_of_PVCS_grid_points[3];
    float relative_coordinates[3];
    float grad_value = 0;
    
    #pragma unroll
    for (uint i=0; i<2; i++)
    {
        #pragma unroll
        for (uint j=0; j<2; j++)
        {
            #pragma unroll
            for (uint k=0; k<2; k++)
            {
                temp_PVCS_index_of_PVCS_grid_points[0] = PVCS_index_of_BEV_grid_points[0] + i;
                temp_PVCS_index_of_PVCS_grid_points[1] = PVCS_index_of_BEV_grid_points[1] + j;
                temp_PVCS_index_of_PVCS_grid_points[2] = PVCS_index_of_BEV_grid_points[2] + k;

                float temp_value = tex3D<float>(PVCSDoseGradTexture, temp_PVCS_index_of_PVCS_grid_points[2], \
                    temp_PVCS_index_of_PVCS_grid_points[1], temp_PVCS_index_of_PVCS_grid_points[0]);
                
                PVCS_coordinates_of_PVCS_grid_points[0] = ((float)(temp_PVCS_index_of_PVCS_grid_points[0]) + half) * voxel_size - phantom_iso0;
                PVCS_coordinates_of_PVCS_grid_points[1] = ((float)(temp_PVCS_index_of_PVCS_grid_points[1]) + half) * voxel_size - phantom_iso1;
                PVCS_coordinates_of_PVCS_grid_points[2] = ((float)(temp_PVCS_index_of_PVCS_grid_points[2]) + half) * voxel_size - phantom_iso2;
                d_PVCS_to_BEV(BEV_coordinates_of_PVCS_grid_points, PVCS_coordinates_of_PVCS_grid_points, zenith, azimuth);

                BEV_coordinates_of_PVCS_grid_points[2] += SAD;
                fluence_map_pixel_size_at_this_depth = fluence_map_pixel_size * BEV_coordinates_of_PVCS_grid_points[2] / SAD;

                BEV_index_of_PVCS_grid_points[0] = middle_point + BEV_coordinates_of_PVCS_grid_points[0] / fluence_map_pixel_size_at_this_depth;
                BEV_index_of_PVCS_grid_points[1] = middle_point + BEV_coordinates_of_PVCS_grid_points[1] / fluence_map_pixel_size_at_this_depth;
                BEV_index_of_PVCS_grid_points[2] = (BEV_coordinates_of_PVCS_grid_points[2] - sampling_start) / sampling_step;

                relative_coordinates[0] = fabsf(BEV_index_of_PVCS_grid_points[0] - idx_x);
                relative_coordinates[1] = fabsf(BEV_index_of_PVCS_grid_points[1] - idx_y);
                relative_coordinates[2] = fabsf(BEV_index_of_PVCS_grid_points[2] - idx_z);
                bool valid = relative_coordinates[0] < 1 && relative_coordinates[1] < 1 && \
                    relative_coordinates[2] < 1;
                float trilinear_interpolation_coefficient = (1 - relative_coordinates[0]) * (1 - relative_coordinates[1]) * \
                    (1 - relative_coordinates[2]);
                grad_value += valid * trilinear_interpolation_coefficient * temp_value;
            }
        }
    }

    // note that d_FCBB_BEV_dose_grad follows logical order: (z, x, y).
    uint result_idx = (idx_z * fluence_map_dimension + idx_x) * fluence_map_dimension + idx_y;
    d_FCBB_BEV_dose_grad[result_idx] = !(idx_x >= fluence_map_dimension - 1 || idx_y >= fluence_map_dimension - 1 \
        || idx_z >= sampling_points) * grad_value;
}


extern "C"
void PVCSDoseBackward(float voxel_size, float phantom_iso[3], \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_pixel_size, uint fluence_map_dimension, \
    float* d_FCBB_BEV_dose_grad, cudaTextureObject_t PVCSDoseGradTexture, \
    cudaStream_t stream)
{
    dim3 blockSize(16, 16, 1); // (x, y, z)
    dim3 gridSize;
    gridSize.x = ceil((float)fluence_map_dimension / blockSize.x);
    gridSize.y = ceil((float)fluence_map_dimension / blockSize.y);
    gridSize.z = ceil((float)sampling_points / blockSize.z);

    d_PVCSDoseBackward<<<gridSize, blockSize, 0, stream>>>(voxel_size, \
        phantom_iso[0], phantom_iso[1], phantom_iso[2], zenith, azimuth, SAD, \
        sampling_start, sampling_end, sampling_points, fluence_map_pixel_size, \
        fluence_map_dimension, d_FCBB_BEV_dose_grad, PVCSDoseGradTexture);
}


__global__ void
d_writeFCBBPVCSDoseGradSurface(cudaSurfaceObject_t surface, float* input, uint dim_x, uint dim_y, uint dim_z)
{
    uint idx_x = blockIdx.x * blockDim.x + threadIdx.x;
    uint idx_y = blockIdx.y * blockDim.y + threadIdx.y;
    uint idx_z = blockIdx.z * blockDim.z + threadIdx.z;

    if(idx_x >= dim_x || idx_y >= dim_y || idx_z >= dim_z)
        return;
    
    uint idx_input = (idx_x * dim_y + idx_y) * dim_z + idx_z;
    float value = input[idx_input];
    surf3Dwrite(value, surface, idx_z*sizeof(float), idx_y, idx_x, cudaBoundaryModeTrap);
}

extern "C"
void writeFCBBPVCSDoseGradSurface(cudaSurfaceObject_t surface, float* input, \
    uint dim_x, uint dim_y, uint dim_z, cudaStream_t stream=0)
{
    dim3 blockSize(1, 8, 8);
    dim3 gridSize(ceil((float)dim_x / blockSize.x), ceil((float)dim_y / blockSize.y), ceil((float)dim_z / blockSize.z));
    d_writeFCBBPVCSDoseGradSurface<<<gridSize, blockSize, 0, stream>>>(surface, input, dim_x, dim_y, dim_z);
}

__global__ void
d_textureMinusCoordinate(cudaTextureObject_t texture, float* d_output, \
    float coord0, float coord1, float coord2)
{
    d_output[0] = tex3D<float>(texture, coord2, coord1, coord0);
}

extern "C"
void textureMinusCoordinate(cudaTextureObject_t texture, float* d_output, \
    float coord0, float coord1, float coord2)
{
    dim3 blockSize(1);
    dim3 gridSize(1);
    d_textureMinusCoordinate<<<gridSize, blockSize>>>(texture, d_output, coord0, coord1, coord2);
}