#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;


void E2E::host_PVCS_to_BEV(float BEV_coords[3], float PVCS_coords[3], float theta, float phi)
{
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    float sin_phi = sin(phi);
    float cos_phi = cos(phi);

    BEV_coords[0] = cos_theta * cos_phi * PVCS_coords[0] + cos_theta * sin_phi * PVCS_coords[1] - sin_theta * PVCS_coords[2];
    BEV_coords[1] = -sin_phi * PVCS_coords[0] + cos_phi * PVCS_coords[1];
    BEV_coords[2] = sin_theta * cos_phi * PVCS_coords[0] + sin_theta * sin_phi * PVCS_coords[1] + cos_theta * PVCS_coords[2];
}

void module_test_PVCS_dose_forward_host(beam& Beam, phantom& Phtm, float* h_PVCS_dose, float* h_BEV_dose)
{
    array<uint, 3> BEV_dimension({Beam.sampling_points, Beam.convolved_fluence_map_dimension[0], \
        Beam.convolved_fluence_map_dimension[1]});
    array<uint, 3> PVCS_dimension({Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch});
    float pixel_size = Beam.pixel_size;
    float sampling_start = Beam.sampling_range[0];
    float sampling_end = Beam.sampling_range[1];
    uint sampling_points = Beam.sampling_points;
    float sampling_step = (sampling_end - sampling_start) / (sampling_points - 1);
    float theta = Beam.zenith;
    float phi = Beam.azimuth;

    float voxel_size = Phtm.voxelSize;
    float PVCS_x_hat[3] = {voxel_size, 0, 0};
    float PVCS_y_hat[3] = {0, voxel_size, 0};
    float PVCS_z_hat[3] = {0, 0, voxel_size};

    float BEV_x_hat[3];
    float BEV_y_hat[3];
    float BEV_z_hat[3];

    host_PVCS_to_BEV(BEV_x_hat, PVCS_x_hat, theta, phi);
    host_PVCS_to_BEV(BEV_y_hat, PVCS_y_hat, theta, phi);
    host_PVCS_to_BEV(BEV_z_hat, PVCS_z_hat, theta, phi);

    float PVCS_origin[3] = {-Beam.isocenter[0], -Beam.isocenter[1], -Beam.isocenter[2]};
    float BEV_origin[3];
    host_PVCS_to_BEV(BEV_origin, PVCS_origin, theta, phi);
    BEV_origin[2] += Beam.SAD;

    constexpr float half = 0.5;
    float x_center = (float)(Beam.convolved_fluence_map_dimension[0]) / 2 - half;
    float y_center = (float)(Beam.convolved_fluence_map_dimension[1]) / 2 - half;
    for (uint i=0; i<PVCS_dimension[0]; i++)
    {
        float BEV_point_i[3];
        BEV_point_i[0] = BEV_origin[0] + (i + half) * BEV_x_hat[0];
        BEV_point_i[1] = BEV_origin[1] + (i + half) * BEV_x_hat[1];
        BEV_point_i[2] = BEV_origin[2] + (i + half) * BEV_x_hat[2];
        uint idx_i = i * PVCS_dimension[1];
        for (uint j=0; j<PVCS_dimension[1]; j++)
        {
            float BEV_point_j[3];
            BEV_point_j[0] = BEV_point_i[0] + (j + half) * BEV_y_hat[0];
            BEV_point_j[1] = BEV_point_i[1] + (j + half) * BEV_y_hat[1];
            BEV_point_j[2] = BEV_point_i[2] + (j + half) * BEV_y_hat[2];
            uint idx_j = (idx_i + j) * PVCS_dimension[2];
            for (uint k=0; k<PVCS_dimension[2]; k++)
            {
                float BEV_point_k[3];
                BEV_point_k[0] = BEV_point_j[0] + (k + half) * BEV_z_hat[0];
                BEV_point_k[1] = BEV_point_j[1] + (k + half) * BEV_z_hat[1];
                BEV_point_k[2] = BEV_point_j[2] + (k + half) * BEV_z_hat[2];
                uint idx_k = idx_j + k;

                float scale = BEV_point_k[2] / Beam.SAD;
                float pixel_size_at_this_depth = pixel_size * scale;

                // normalize
                float normalized_BEV_coords[3] = {(BEV_point_k[2] - sampling_start) / sampling_step + half, \
                    BEV_point_k[0] / pixel_size_at_this_depth + x_center, \
                    BEV_point_k[1] / pixel_size_at_this_depth + y_center};

                h_PVCS_dose[idx_k] = read_phantom(normalized_BEV_coords, BEV_dimension, h_BEV_dose);
            }
        }
    }
}


extern "C"
void writeSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::module_test_PVCS_dose_forward(beam& Beam, phantom& Phtm)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolved_fluence_map_size = convolved_fluence_map_dimension * convolved_fluence_map_dimension;
    uint sampling_points = Beam.sampling_points;
    uint FCBB_BEV_dose_size = convolved_fluence_map_size * sampling_points;

    array<uint, 3> phantom_dimension{Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    uint FCBB_PVCS_dose_size = phantom_dimension[0] * phantom_dimension[1] * phantom_dimension[2];

    float* h_FCBB_BEV_dose = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));
    float* h_FCBB_PVCS_dose = (float*)malloc(FCBB_PVCS_dose_size*sizeof(float));
    float* h_FCBB_PVCS_dose_d = (float*)malloc(FCBB_PVCS_dose_size*sizeof(float));
    // initialize
    uint norm = 65536;
    for (uint i=0; i<FCBB_BEV_dose_size; i++)
        h_FCBB_BEV_dose[i] = (float)(rand() % norm) / (norm - 1);
    
    // // for debug purposes
    // for (uint i=0; i<sampling_points; i++)
    // {
    //     uint idx_i = i * convolved_fluence_map_size;
    //     for (uint j=0; j<convolved_fluence_map_size; j++)
    //     {
    //         uint idx = idx_i + j;
    //         h_FCBB_BEV_dose[idx] = i;
    //     }
    // }
    
    // load to device
    float* d_FCBB_BEV_dose;
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_BEV_dose, FCBB_BEV_dose_size*sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_FCBB_BEV_dose, h_FCBB_BEV_dose, \
        FCBB_BEV_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    dim3 blockSize(1, 16, 16);
    dim3 gridSize(sampling_points / blockSize.x, convolved_fluence_map_dimension / blockSize.y, \
        convolved_fluence_map_dimension / blockSize.z);
    writeSurface(gridSize, blockSize, Beam.FCBB_BEV_dose_texture, d_FCBB_BEV_dose);
    checkCudaErrors(cudaFree(d_FCBB_BEV_dose));

    // device compute
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "device compute begin" << endl;

    cudaEventRecord(start);
    Beam.PVCS_dose_forward(Phtm, 0);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << "ms" << endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // offload to cpu
    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose_d, Beam.d_FCBB_PVCS_dose, \
        FCBB_PVCS_dose_size*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    cout << "Host compute" << endl;
    module_test_PVCS_dose_forward_host(Beam, Phtm, h_FCBB_PVCS_dose, h_FCBB_BEV_dose);

    // compare
    float mse = 0;
    for (uint i=0; i<FCBB_PVCS_dose_size; i++)
    {
        float diff = h_FCBB_PVCS_dose[i] - h_FCBB_PVCS_dose_d[i];
        mse += diff * diff;
    }
    mse /= FCBB_PVCS_dose_size;
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<FCBB_PVCS_dose_size; i++)
        base += h_FCBB_PVCS_dose[i] * h_FCBB_PVCS_dose[i];
    base /= FCBB_PVCS_dose_size;
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // // for debug purposes
    // string MT_PVCS_host{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/MT_PVCS_host.dat"};
    // string MT_PVCS_device{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/MT_PVCS_device.dat"};
    // ofstream outFile(MT_PVCS_host);
    // outFile.write((char*)h_FCBB_PVCS_dose, FCBB_PVCS_dose_size*sizeof(float));
    // outFile.close();
    // outFile.open(MT_PVCS_device);
    // outFile.write((char*)h_FCBB_PVCS_dose_d, FCBB_PVCS_dose_size*sizeof(float));
    // outFile.close();

    // cleanup
    free(h_FCBB_BEV_dose);
    free(h_FCBB_PVCS_dose);
}