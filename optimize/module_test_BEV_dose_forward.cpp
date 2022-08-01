#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <assert.h>
#include "geom.h"
#include "args.h"
#include "optim.h"

using namespace E2E;
using namespace std;

float read_phantom(float coords[3], array<uint, 3>& phantom_dimension, float* phantom)
{
    // coords are unitless, the same as phantom_dimension
    constexpr float shift = 0.5;
    float coords_normalized[3] = {coords[0] - shift, coords[1] - shift, coords[2] - shift}; // shift 0.5 to 0
    int coords_base[3] = {floor(coords_normalized[0]), floor(coords_normalized[1]), floor(coords_normalized[2])};
    float coords_shift[3] = {coords_normalized[0] - coords_base[0], \
        coords_normalized[1] - coords_base[1], coords_normalized[2] - coords_base[2]};
    
    float value = 0;
    for (uint i=0; i<2; i++)
    {
        int coords_i = coords_base[0] + i;
        if (coords_i < 0 || coords_i >= phantom_dimension[0]) continue;
        float coeff_i = abs(1 - i - coords_shift[0]);
        for (uint j=0; j<2; j++)
        {
            int coords_j = coords_base[1] + j;
            if (coords_j < 0 || coords_j >= phantom_dimension[1]) continue;
            coords_j = coords_i * phantom_dimension[1] + coords_j;
            float coeff_j = abs(1 - j - coords_shift[1]) * coeff_i;
            for (uint k=0; k<2; k++)
            {
                int coords_k = coords_base[2] + k;
                if (coords_k < 0 || coords_k >= phantom_dimension[2]) continue;
                coords_k = coords_j * phantom_dimension[2] + coords_k;
                float coeff_k = abs(1 - k - coords_shift[2]) * coeff_j;
                value += coeff_k * phantom[coords_k];
            }
        }
    }
    return value;
}

float read_depth_dose(float* h_depth_dose, float coord, uint size)
{
    // h_depth_dose stores the depth dose data
    // size is the number of elements of h_depth_dose
    float coord_norm = coord - 0.5;
    int coord_base = floor(coord_norm);
    float coord_shift = coord_norm - coord_base;
    float value = 0;
    for (uint i=0; i<2; i++)
    {
        int coord_i = coord_base + i;
        if (coord_i < 0 || coord_i >= size)
            continue;
        value += h_depth_dose[coord_i] * abs(1. - i - coord_shift);
    }
    return value;
}

extern "C"
void moduleTestHostTrilinear(cudaTextureObject_t tex, float* d_result, uint dim0, uint dim1, \
    float level2, float shift0, float shift1);

void E2E::module_test_host_triliner(phantom& Phtm)
{
    uint dim0 = 111;
    uint dim1 = 111;
    float level2 = 0.62;
    float shift0 = 0.16;
    float shift1 = 0.28;
    uint size = dim0 * dim1;

    float* h_result = (float*)malloc(size*sizeof(float));
    float* h_result_d = (float*)malloc(size*sizeof(float));
    float* d_result;
    checkCudaErrors(cudaMalloc((void**)&d_result, size*sizeof(float)));

    // device compute
    moduleTestHostTrilinear(Phtm.tex, d_result, dim0, dim1, level2, shift0, shift1);

    // offload
    checkCudaErrors(cudaMemcpy(h_result_d, d_result, size*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    assert(Phtm.pitchPadding);
    float* h_HU = Phtm.h_HU;
    array<uint, 3> phantom_dimension({Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch});
    float coords[3];
    coords[2] = phantom_dimension[2] * level2;
    for (uint i=0; i<dim0; i++)
    {
        for (uint j=0; j<dim1; j++)
        {
            uint result_idx = i * dim1 + j;
            coords[0] = ((float)i / dim0 + shift0) * phantom_dimension[0];
            coords[1] = ((float)j / dim1 + shift1) * phantom_dimension[1];
            h_result[result_idx] = read_phantom(coords, phantom_dimension, h_HU);
        }
    }

    // compare
    float MSE = 0;
    for (uint i=0; i<size; i++)
    {
        float diff = h_result_d[i] - h_result[i];
        MSE += diff * diff;
    }
    MSE /= size;
    MSE = sqrt(MSE);

    float base = 0;
    for (uint i=0; i<size; i++)
        base += h_result[i] * h_result[i];
    base /= size;
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << MSE << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_result);
    free(h_result_d);
    checkCudaErrors(cudaFree(d_result));
}

extern "C"
void moduleTestHostLinear(cudaTextureObject_t tex, float* d_result, uint dim, float shift);

void E2E::module_test_host_linear()
{
    uint n_samples = 107;
    uint num_depths = (*FCBB6MeV).num_depths;
    float shift = 0.19;
    float* h_result = (float*)malloc(n_samples*sizeof(float));
    float* h_result_d = (float*)malloc(n_samples*sizeof(float));
    float* d_result;
    checkCudaErrors(cudaMalloc((void**)&d_result, n_samples*sizeof(float)));

    // device compute
    moduleTestHostLinear((*FCBB6MeV).tex, d_result, n_samples, shift);

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_result_d, d_result, n_samples*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    float* h_doses = (*FCBB6MeV).doses;
    for (uint i=0; i<n_samples; i++)
    {
        float coord = ((float)i / n_samples + shift) * num_depths;
        h_result[i] = read_depth_dose(h_doses, coord, num_depths);
    }

    // compare
    float mse = 0;
    for (uint i=0; i<n_samples; i++)
    {
        float diff = h_result_d[i] - h_result[i];
        mse += diff * diff;
    }
    mse /= n_samples;
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<n_samples; i++)
        base += h_result[i] * h_result[i];
    base /= n_samples;
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_result);
    free(h_result_d);
    checkCudaErrors(cudaFree(d_result));
}

void host_BEV_to_PVCS(float PVCS_coords[3], float BEV_coords[3], float theta, float phi)
{
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    float sin_phi = sin(phi);
    float cos_phi = cos(phi);

    PVCS_coords[0] = cos_theta * cos_phi * BEV_coords[0] - sin_phi * BEV_coords[1] + sin_theta * cos_phi * BEV_coords[2];
    PVCS_coords[1] = cos_theta * sin_phi * BEV_coords[0] + cos_phi * BEV_coords[1] + sin_theta * sin_phi * BEV_coords[2];
    PVCS_coords[2] = - sin_theta * BEV_coords[0] + cos_theta * BEV_coords[2];
}

void module_test_BEV_dose_forward_host(beam& Beam, phantom& Phtm, float* h_FCBB_BEV_dose, float* convolved_fluence_map)
{
    // h_FCBB_BEV_dose is used for storing results
    assert(Phtm.pitchPadding);
    array<uint, 3> phantom_dimension({Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch});
    array<float, 3> phantom_isocenter({Beam.isocenter[0], Beam.isocenter[1], Beam.isocenter[2]});
    float phantom_voxel_size = Phtm.voxelSize;
    float* h_HU = Phtm.h_HU;

    float* depthDose = (*FCBB6MeV).doses;
    
    uint depthDosePoints = 400;
    uint sampling_points = Beam.sampling_points;
    float max_depth = (*FCBB6MeV).depths[depthDosePoints-1];
    float fluence_map_pixel_size = Beam.pixel_size;
    uint convolved_fluence_map_dimension = Beam.convolved_fluence_map_dimension[0];
    float sampling_range_start = Beam.sampling_range[0];
    float sampling_range_end = Beam.sampling_range[1];

    float theta = Beam.zenith;
    float phi = Beam.azimuth;
    float BEV_source_coords[3] = {0, 0, -Beam.SAD};
    float PVCS_source_coords[3];
    host_BEV_to_PVCS(PVCS_source_coords, BEV_source_coords, theta, phi);
    PVCS_source_coords[0] += phantom_isocenter[0];
    PVCS_source_coords[1] += phantom_isocenter[1];
    PVCS_source_coords[2] += phantom_isocenter[2];

    float BEV_begin[3];
    float PVCS_begin[3];
    float PVCS_coords[3];
    float PVCS_coords_normalized[3];
    float x_center = (float)(convolved_fluence_map_dimension - 2) / 2;
    float y_center = (float)(convolved_fluence_map_dimension - 2) / 2;
    float sampling_step = (sampling_range_end - sampling_range_start) / sampling_points;

    for (uint i=0; i<convolved_fluence_map_dimension; i++)
    {
        for (uint j=0; j<convolved_fluence_map_dimension; j++)
        {
            BEV_begin[0] = ((float)i - x_center) * fluence_map_pixel_size;
            BEV_begin[1] = ((float)j - y_center) * fluence_map_pixel_size;
            BEV_begin[2] = Beam.SAD;
            host_BEV_to_PVCS(PVCS_begin, BEV_begin, theta, phi);
            float radiological_path_step = sqrt(BEV_begin[0] * BEV_begin[0] + \
                BEV_begin[1] * BEV_begin[1] + BEV_begin[2] * BEV_begin[2]) * sampling_step / Beam.SAD;
            
            float distance = 0;
            float scale = 0;
            float radiological_path_length = 0;

            uint idx_ij = i * convolved_fluence_map_dimension + j;
            for (uint k=0; k<sampling_points; k++)
            {
                distance = sampling_range_start + k * sampling_step;
                scale = distance / Beam.SAD;
                PVCS_coords[0] = PVCS_source_coords[0] + PVCS_begin[0] * scale;
                PVCS_coords[1] = PVCS_source_coords[1] + PVCS_begin[1] * scale;
                PVCS_coords[2] = PVCS_source_coords[2] + PVCS_begin[2] * scale;

                PVCS_coords_normalized[0] = PVCS_coords[0] / phantom_voxel_size;
                PVCS_coords_normalized[1] = PVCS_coords[1] / phantom_voxel_size;
                PVCS_coords_normalized[2] = PVCS_coords[2] / phantom_voxel_size;

                float HU = read_phantom(PVCS_coords_normalized, phantom_dimension, h_HU);

                radiological_path_length += HU * radiological_path_step;
                float normalized_radiological_path_length = radiological_path_length / max_depth * depthDosePoints;
                float dose = read_depth_dose(depthDose, normalized_radiological_path_length, depthDosePoints);

                uint idx = k * convolved_fluence_map_dimension * convolved_fluence_map_dimension + idx_ij;
                h_FCBB_BEV_dose[idx] = dose * convolved_fluence_map[idx_ij] / (scale * scale);
            }
        }
    }
}

extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::module_test_BEV_dose_forward(beam& Beam, phantom& Phtm)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolved_fluence_map_size = convolved_fluence_map_dimension * convolved_fluence_map_dimension;

    // convolved_fluence_map initialization
    float* h_convolved_fluence_map = (float*)malloc(convolved_fluence_map_size*sizeof(float));
    for (uint i=0; i<convolved_fluence_map_size; i++)
        h_convolved_fluence_map[i] = 1;
    checkCudaErrors(cudaMemcpy(Beam.d_convolved_fluence_map, h_convolved_fluence_map, \
        convolved_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));

    uint sampling_points = Beam.sampling_points;
    uint FCBB_BEV_dose_size = convolved_fluence_map_size * sampling_points;
    float* h_FCBB_BEV_dose = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));
    float* h_FCBB_BEV_dose_d = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));
    float* d_FCBB_BEV_dose;
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_BEV_dose, FCBB_BEV_dose_size*sizeof(float)));

    // device compute
    Beam.BEV_dose_forward(Phtm, FCBB6MeV, 0);

    // offload to CPU
    // FCBB_BEV_dose_texture is ordered in (z, x, y). cudaExtent order: (y, x, z)
    dim3 blockSize(16, 16, 1);
    dim3 gridSize(convolved_fluence_map_dimension / blockSize.x, \
        convolved_fluence_map_dimension / blockSize.y, sampling_points / blockSize.z);
    readSurface(gridSize, blockSize, Beam.FCBB_BEV_dose_surface, d_FCBB_BEV_dose);
    checkCudaErrors(cudaMemcpy(h_FCBB_BEV_dose_d, d_FCBB_BEV_dose, \
        FCBB_BEV_dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(d_FCBB_BEV_dose));

    // host compute
    cout << "Host compute" << endl;
    module_test_BEV_dose_forward_host(Beam, Phtm, h_FCBB_BEV_dose, h_convolved_fluence_map);

    // compare
    float mse = 0;
    for (uint i=0; i<FCBB_BEV_dose_size; i++)
    {
        float diff = h_FCBB_BEV_dose_d[i] - h_FCBB_BEV_dose[i];
        mse += diff * diff;
    }
    mse /= FCBB_BEV_dose_size;
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<FCBB_BEV_dose_size; i++)
        base += h_FCBB_BEV_dose[i] * h_FCBB_BEV_dose[i];
    base /= FCBB_BEV_dose_size;
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // string host_result_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/BEV_dose_host.dat"};
    // string device_result_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/BEV_dose_device.dat"};
    // ofstream outFile(host_result_path);
    // outFile.write((char*)h_FCBB_BEV_dose, FCBB_BEV_dose_size*sizeof(float));
    // outFile.close();
    // outFile.open(device_result_path);
    // outFile.write((char*)h_FCBB_BEV_dose_d, FCBB_BEV_dose_size*sizeof(float));
    // outFile.close();

    // cleanup
    free(h_convolved_fluence_map);
    free(h_FCBB_BEV_dose);
    free(h_FCBB_BEV_dose_d);
}