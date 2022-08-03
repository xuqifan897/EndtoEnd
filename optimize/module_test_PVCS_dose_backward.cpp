#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

// void E2E::write_phantom(float coords[3], std::array<uint, 3>& phantom_dimension, float* phantom, float value)
// {
//     constexpr float shift = 0.5;
//     float coords_normalized[3] = {coords[0] - shift, coords[1] - shift, coords[2] - shift}; // shift 0.5 to 0
//     int coords_base[3] = {floor(coords_normalized[0]), floor(coords_normalized[1]), floor(coords_normalized[2])};
//     float coords_shift[3] = {coords_normalized[0] - coords_base[0], \
//         coords_normalized[1] - coords_base[1], coords_normalized[2] - coords_base[2]};
    
//     for (uint i=0; i<2; i++)
//     {
//         int coords_i = coords_base[0] + i;
//         if (coords_i < 0 || coords_i >= phantom_dimension[0]) continue;
//         float coeff_i = abs(1 - i - coords_shift[0]);
//         for (uint j=0; j<2; j++)
//         {
//             int coords_j = coords_base[1] + j;
//             if (coords_j < 0 || coords_j >= phantom_dimension[1]) continue;
//             coords_j = coords_i * phantom_dimension[1] + coords_j;
//             float coeff_j = abs(1 - j - coords_shift[1]) * coeff_i;
//             for (uint k=0; k<2; k++)
//             {
//                 int coords_k = coords_base[2] + k;
//                 if (coords_k < 0 || coords_k >= phantom_dimension[2]) continue;
//                 coords_k = coords_j * phantom_dimension[2] + coords_k;
//                 float coeff_k = abs(1 - k - coords_shift[2]) * coeff_j;
//                 phantom[coords_k] += coeff_k * value;
//             }
//         }
//     }
// }

void E2E::write_phantom(float coords[3], array<uint, 3>& phantom_dimension, float* phantom, float value)
{
    constexpr float half = 0.5;
    float normalized_coords[3] = {coords[0] - half, coords[1] - half, coords[2] - half};
    int coords_base[3] = {floor(normalized_coords[0]), floor(normalized_coords[1]), floor(normalized_coords[2])};
    float coords_shift[3] = {normalized_coords[0] - coords_base[0], normalized_coords[1] - coords_base[1], \
        normalized_coords[2] - coords_base[2]};
    
    for (uint i=0; i<2; i++)
    {
        int coords_i = coords_base[0] + i;
        if (coords_i < 0 || coords_i >= phantom_dimension[0]) continue;
        coords_i *= phantom_dimension[1];
        float coeff_i = abs(1 - i - coords_shift[0]);

        for (uint j=0; j<2; j++)
        {
            int coords_j = coords_base[1] + j;
            if (coords_j < 0 || coords_j >= phantom_dimension[1]) continue;
            coords_j = (coords_i + coords_j) * phantom_dimension[2];
            float coeff_j = coeff_i * abs(1 - j - coords_shift[1]);

            for (uint k=0; k<2; k++)
            {
                int coords_k = coords_base[2] + k;
                if (coords_k < 0 || coords_k >= phantom_dimension[2]) continue;
                coords_k += coords_j;
                float coeff_k = coeff_j * abs(1 - k - coords_shift[2]);
                phantom[coords_k] += coeff_k * value;
            }
        }
    }
}

void module_test_PVCS_dose_backward_host(beam& Beam, phantom& Phtm, float* h_FCBB_BEV_dose_grad, float* h_FCBB_PVCS_dose_grad)
{
    float theta = Beam.zenith;
    float phi = Beam.azimuth;
    array<float, 3> isocenter = Beam.isocenter;
    float pixel_size = Beam.pixel_size;
    float voxel_size = Phtm.voxelSize;
    float SAD = Beam.SAD;
    float sampling_step = (Beam.sampling_range[1] - Beam.sampling_range[0]) / (Beam.sampling_points - 1);

    float phantom_dimension[3] = {Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    array<uint, 3> BEV_dose_dimension({Beam.sampling_points, Beam.convolved_fluence_map_dimension[0], \
        Beam.convolved_fluence_map_dimension[1]});
    float PVCS_origin[3] = {-isocenter[0], -isocenter[1], -isocenter[2]};
    float BEV_origin[3];
    host_PVCS_to_BEV(BEV_origin, PVCS_origin, theta, phi);
    BEV_origin[2] += SAD;

    float PVCS_x_hat[3] = {voxel_size, 0, 0};
    float PVCS_y_hat[3] = {0, voxel_size, 0};
    float PVCS_z_hat[3] = {0, 0, voxel_size};

    float BEV_x_hat[3];
    float BEV_y_hat[3];
    float BEV_z_hat[3];

    host_PVCS_to_BEV(BEV_x_hat, PVCS_x_hat, theta, phi);
    host_PVCS_to_BEV(BEV_y_hat, PVCS_y_hat, theta, phi);
    host_PVCS_to_BEV(BEV_z_hat, PVCS_z_hat, theta, phi);

    constexpr float half = 0.5;
    float x_center = (float)(Beam.convolved_fluence_map_dimension[0]) / 2 - half;
    float y_center = (float)(Beam.convolved_fluence_map_dimension[1]) / 2 - half;
    for (uint i=0; i<phantom_dimension[0]; i++)
    {
        float BEV_coords_i[3];
        BEV_coords_i[0] = BEV_origin[0] + (i + half) * BEV_x_hat[0];
        BEV_coords_i[1] = BEV_origin[1] + (i + half) * BEV_x_hat[1];
        BEV_coords_i[2] = BEV_origin[2] + (i + half) * BEV_x_hat[2];
        uint idx_i = i * phantom_dimension[1];

        for (uint j=0; j<phantom_dimension[1]; j++)
        {
            float BEV_coords_j[3];
            BEV_coords_j[0] = BEV_coords_i[0] + (j + half) * BEV_y_hat[0];
            BEV_coords_j[1] = BEV_coords_i[1] + (j + half) * BEV_y_hat[1];
            BEV_coords_j[2] = BEV_coords_i[2] + (j + half) * BEV_y_hat[2];
            uint idx_j = (idx_i + j) * phantom_dimension[2];

            for (uint k=0; k<phantom_dimension[2]; k++)
            {
                float BEV_coords_k[3];
                BEV_coords_k[0] = BEV_coords_j[0] + (k + half) * BEV_z_hat[0];
                BEV_coords_k[1] = BEV_coords_j[1] + (k + half) * BEV_z_hat[1];
                BEV_coords_k[2] = BEV_coords_j[2] + (k + half) * BEV_z_hat[2];
                uint idx_k = idx_j + k;

                float scale = BEV_coords_k[2] / SAD;
                float pixel_size_at_this_depth = pixel_size * scale;
                float normalized_BEV_coords[3] = {(BEV_coords_k[2] - Beam.sampling_range[0]) / sampling_step + half, \
                    BEV_coords_k[0] / pixel_size_at_this_depth + x_center, \
                    BEV_coords_k[1] / pixel_size_at_this_depth + y_center};
                write_phantom(normalized_BEV_coords, BEV_dose_dimension, h_FCBB_BEV_dose_grad, h_FCBB_PVCS_dose_grad[idx_k]);
            }
        }
    }
}

extern "C"
void writeSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::module_test_PVCS_dose_backward(beam& Beam, phantom& Phtm)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolved_fluence_map_size = convolved_fluence_map_dimension * convolved_fluence_map_dimension;
    uint sampling_points = Beam.sampling_points;
    uint FCBB_BEV_dose_size = convolved_fluence_map_size * sampling_points;

    array<uint, 3> phantom_dimension{Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    uint FCBB_PVCS_dose_size = phantom_dimension[0] * phantom_dimension[1] * phantom_dimension[2];

    float* h_FCBB_BEV_dose_grad = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));
    float* h_FCBB_BEV_dose_grad_d = (float*)malloc(FCBB_BEV_dose_size*sizeof(float));

    float* h_FCBB_PVCS_dose_grad = (float*)malloc(FCBB_PVCS_dose_size*sizeof(float));
    // initialize
    uint norm = 65536;
    for (uint i=0; i<FCBB_PVCS_dose_size; i++)
        h_FCBB_PVCS_dose_grad[i] = (float)(rand() % norm) / (norm - 1);
    
    // load to device
    float* d_FCBB_PVCS_dose_grad;
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_PVCS_dose_grad, FCBB_PVCS_dose_size*sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_FCBB_PVCS_dose_grad, h_FCBB_PVCS_dose_grad, \
        FCBB_PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(phantom_dimension[0] / blockSize.x, phantom_dimension[1] / blockSize.y, \
        phantom_dimension[2] / blockSize.z);
    writeSurface(gridSize, blockSize, beam::FCBB_PVCS_dose_grad_surface, d_FCBB_PVCS_dose_grad);
    checkCudaErrors(cudaFree(d_FCBB_PVCS_dose_grad));

    // device compute
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "Device compute" << endl;
    cudaEventRecord(start);
    Beam.PVCS_dose_backward(Phtm, 0);
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << endl;

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_FCBB_BEV_dose_grad_d, Beam.d_FCBB_BEV_dose_grad, \
        FCBB_BEV_dose_size*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    cout << "Host compute" << endl;
    for (uint i=0; i<FCBB_BEV_dose_size; i++)
        h_FCBB_BEV_dose_grad[i] = 0;
    module_test_PVCS_dose_backward_host(Beam, Phtm, h_FCBB_BEV_dose_grad, h_FCBB_PVCS_dose_grad);

    // compare
    float mse = 0;
    float base = 0;
    for (uint i=0; i<sampling_points; i++)
    {
        uint idx_i = i * convolved_fluence_map_dimension;
        for (uint j=0; j<convolved_fluence_map_dimension-1; j++)
        {
            uint idx_j = (idx_i + j) * convolved_fluence_map_dimension;
            for (uint k=0; k<convolved_fluence_map_dimension-1; k++)
            {
                uint idx_k = idx_j + k;
                float diff = h_FCBB_BEV_dose_grad[idx_k] - h_FCBB_BEV_dose_grad_d[idx_k];
                mse += diff * diff;
                base += h_FCBB_BEV_dose_grad[idx_k] * h_FCBB_BEV_dose_grad[idx_k];
            }
        }
    }
    mse /= FCBB_BEV_dose_size;
    base /= FCBB_BEV_dose_size;
    mse = sqrt(mse);
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // // for debug purposes
    // string host_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/BEV_grad_host.dat"};
    // string device_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_debug/BEV_grad_device.dat"};
    // ofstream outFile(host_path);
    // outFile.write((char*)h_FCBB_BEV_dose_grad, FCBB_BEV_dose_size*sizeof(float));
    // outFile.close();
    // outFile.open(device_path);
    // outFile.write((char*)h_FCBB_BEV_dose_grad_d, FCBB_BEV_dose_size*sizeof(float));
    // outFile.close();

    // cleanup
    free(h_FCBB_BEV_dose_grad);
    free(h_FCBB_BEV_dose_grad_d);
    free(h_FCBB_PVCS_dose_grad);
}