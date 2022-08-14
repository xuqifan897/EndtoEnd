#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

void module_test_BEV_dose_backward_host(beam& Beam, phantom& Phtm, float* h_convolved_fluence_grad, float* h_BEV_dose_grad)
{
    // h_convolved_fluence_grad is used to store the host result
    assert(Phtm.pitchPadding);
    array<uint, 3> phantom_dimension({Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch});
    array<float, 3> phantom_isocenter({Beam.isocenter[0], Beam.isocenter[1], Beam.isocenter[2]});
    float phantom_voxel_size = Phtm.voxelSize;
    float* h_HU = Phtm.h_HU;

    float* depthDose = (*FCBB6MeV).doses;
    
    uint depthDosePoints = (*FCBB6MeV).num_depths;
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
    float sampling_step = (sampling_range_end - sampling_range_start) / (sampling_points - 1);

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
            h_convolved_fluence_grad[idx_ij] = 0;
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
                h_convolved_fluence_grad[idx_ij] += dose * h_BEV_dose_grad[idx] / (scale * scale);
            }
        }
    }
}

void E2E::module_test_BEV_dose_backward(beam& Beam, phantom& Phtm)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolved_fluence_map_size = convolved_fluence_map_dimension * convolved_fluence_map_dimension;

    float* h_convolved_fluence_map_grad = (float*)malloc(convolved_fluence_map_size*sizeof(float));
    float* h_convolved_fluence_map_grad_d = (float*)malloc(convolved_fluence_map_size*sizeof(float));

    uint sampling_points = Beam.sampling_points;
    uint FCBB_BEV_dose_size = convolved_fluence_map_size * sampling_points;

    float* h_FCBB_BEV_dose_grad = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));
    // initialize
    uint norm = 65536;
    for (uint i=0; i<FCBB_BEV_dose_size; i++)
        h_FCBB_BEV_dose_grad[i] = (float)(rand() % norm) / (norm - 1);
    checkCudaErrors(cudaMemcpy(Beam.d_FCBB_BEV_dose_grad, h_FCBB_BEV_dose_grad, \
        FCBB_BEV_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    
    // device compute
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "device compute begin" << endl;

    cudaEventRecord(start);
    Beam.BEV_dose_backward(Phtm, FCBB6MeV, 0);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << "ms" << endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_convolved_fluence_map_grad_d, Beam.d_convolved_fluence_map_grad, \
        convolved_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
    
    // host compute
    cout << "Host compute" << endl;
    module_test_BEV_dose_backward_host(Beam, Phtm, h_convolved_fluence_map_grad, h_FCBB_BEV_dose_grad);

    // compare
    float mse = 0;
    for (uint i=0; i<convolved_fluence_map_dimension-1; i++)
    {
        for (uint j=0; j<convolved_fluence_map_dimension-1; j++)
        {
            uint idx = i * convolved_fluence_map_dimension + j;
            float diff = h_convolved_fluence_map_grad[idx] - h_convolved_fluence_map_grad_d[idx];
            mse += diff * diff;
        }
    }
    mse /= (convolved_fluence_map_dimension - 1) * (convolved_fluence_map_dimension - 1);
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<convolved_fluence_map_dimension-1; i++)
    {
        for (uint j=0; j<convolved_fluence_map_dimension-1; j++)
        {
            uint idx = i * convolved_fluence_map_dimension + j;
            base += h_convolved_fluence_map_grad[idx] * h_convolved_fluence_map_grad[idx];
        }
    }
    base /= (convolved_fluence_map_dimension - 1) * (convolved_fluence_map_dimension-1);
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // string CF_grad_host_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/CF_grad_host.dat"};
    // string CF_grad_device_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/CF_grad_device.dat"};
    // ofstream outFile(CF_grad_host_path);
    // outFile.write((char*)h_convolved_fluence_map_grad, convolved_fluence_map_size*sizeof(float));
    // outFile.close();
    // outFile.open(CF_grad_device_path);
    // outFile.write((char*)h_convolved_fluence_map_grad_d, convolved_fluence_map_size*sizeof(float));
    // outFile.close();

    // cleanup
    free(h_convolved_fluence_map_grad);
    free(h_convolved_fluence_map_grad_d);
    free(h_FCBB_BEV_dose_grad);
}