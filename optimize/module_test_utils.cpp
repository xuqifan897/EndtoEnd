#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

void E2E::module_test_dose_sum(std::vector<beam>& beams, phantom& Phtm)
{
    uint PVCS_dose_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_PVCS_dose_temp = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_PVCS_dose_result = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_PVCS_dose_result_d = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* d_PVCS_dose_result;
    checkCudaErrors(cudaMalloc((void**)&d_PVCS_dose_result, PVCS_dose_size*sizeof(float)));

    // host result initialization
    for (uint i=0; i<PVCS_dose_size; i++)
        h_PVCS_dose_result[i] = 0;

    uint norm = 65535;
    for (uint i=0; i<beams.size(); i++)
    {
        // initialize h_PVCS_dose_temp
        for (uint j=0; j<PVCS_dose_size; j++)
            h_PVCS_dose_temp[j] = (float)(rand() % norm) / (norm - 1);
        
        // copy to device
        checkCudaErrors(cudaMemcpy(beams[i].d_FCBB_PVCS_dose, h_PVCS_dose_temp, \
            PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
        
        // host compute
        for (uint j=0; j<PVCS_dose_size; j++)
            h_PVCS_dose_result[j] += h_PVCS_dose_temp[j];
    }

    // device compute
    float** h_sources = (float**)malloc(beams.size()*sizeof(float*));
    for (uint i=0; i<beams.size(); i++)
        h_sources[i] = beams[i].d_FCBB_PVCS_dose;
    float** d_sources;
    checkCudaErrors(cudaMalloc((void***)&d_sources, beams.size()*sizeof(float*)));
    checkCudaErrors(cudaMemcpy(d_sources, h_sources, beams.size()*sizeof(float*), cudaMemcpyHostToDevice));

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "Device compute" << endl;
    cudaEventRecord(start);
    dose_sum(beams, Phtm, &d_PVCS_dose_result, d_sources, 0);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << "ms" << endl;

    // load to CPU
    checkCudaErrors(cudaMemcpy(h_PVCS_dose_result_d, d_PVCS_dose_result, \
        PVCS_dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    
    // compare
    float mse = 0;
    for (uint i=0; i<PVCS_dose_size; i++)
    {
        float diff = h_PVCS_dose_result[i] - h_PVCS_dose_result_d[i];
        mse += diff * diff;
    }
    mse /= PVCS_dose_size;
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<PVCS_dose_size; i++)
        base += h_PVCS_dose_result[i] * h_PVCS_dose_result[i];
    base /= PVCS_dose_size;
    base = sqrt(base);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_PVCS_dose_temp);
    free(h_PVCS_dose_result);
    free(h_PVCS_dose_result_d);
    checkCudaErrors(cudaFree(d_PVCS_dose_result));
    free(h_sources);
    checkCudaErrors(cudaFree(d_sources));
}

extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::module_test_PVCS_dose_grad(phantom& Phtm)
{
    uint PVCS_dose_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_element_wise_loss = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_PVCS_dose_grad = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_PVCS_total_dose = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_element_wise_loss_d = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* h_PVCS_dose_grad_d = (float*)malloc(PVCS_dose_size*sizeof(float));
    float* d_element_wise_loss;
    float* d_PVCS_dose_grad;
    float* d_PVCS_total_dose;
    checkCudaErrors(cudaMalloc((void**)&d_element_wise_loss, PVCS_dose_size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_PVCS_dose_grad, PVCS_dose_size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_PVCS_total_dose, PVCS_dose_size*sizeof(float)));

    // initialize
    uint norm = 65536;
    for (uint i=0; i<PVCS_dose_size; i++)
    {
        h_PVCS_total_dose[i] = (float)(rand() % norm) / (norm - 1);
        Phtm.h_PTVweight[i] = (float)(rand() % norm) / (norm - 1);
        Phtm.h_PTVtarget[i] = (float)(rand() % norm) / (norm - 1);
        Phtm.h_OARweight[i] = (float)(rand() % norm) / (norm - 1);
        Phtm.h_OARtarget[i] = (float)(rand() % norm) / (norm - 1);
    }

    // load to device
    checkCudaErrors(cudaMemcpy(d_PVCS_total_dose, h_PVCS_total_dose, \
        PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVweight, Phtm.h_PTVweight, \
        PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVtarget, Phtm.h_PTVtarget, \
        PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARweight, Phtm.h_OARweight, \
        PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARtarget, Phtm.h_OARtarget, \
        PVCS_dose_size*sizeof(float), cudaMemcpyHostToDevice));

    // device compute
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "Device compute" << endl;
    cudaEventRecord(start);
    beam::calc_FCBB_PVCS_dose_grad(Phtm, &d_element_wise_loss, d_PVCS_total_dose, 0);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << "ms" << endl;

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_element_wise_loss_d, d_element_wise_loss, \
        PVCS_dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(Phtm.dimension[0] / blockSize.x, Phtm.dimension[1] / blockSize.y, Phtm.pitch / blockSize.z);
    readSurface(gridSize, blockSize, beam::FCBB_PVCS_dose_grad_surface, d_PVCS_dose_grad);
    checkCudaErrors(cudaMemcpy(h_PVCS_dose_grad_d, d_PVCS_dose_grad, PVCS_dose_size*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    cout << "Host compute" << endl;
    for (uint i=0; i<PVCS_dose_size; i++)
    {   
        float PTV_diff = h_PVCS_total_dose[i] - Phtm.h_PTVtarget[i];
        float OAR_diff = h_PVCS_total_dose[i] - Phtm.h_OARtarget[i];
        bool OAR_diff_flag = OAR_diff > 0;
        h_PVCS_dose_grad[i] = 2 * Phtm.h_PTVweight[i] * PTV_diff + \
            2 * Phtm.h_OARweight[i] * OAR_diff * OAR_diff_flag;
        h_element_wise_loss[i] = Phtm.h_PTVweight[i] * PTV_diff * PTV_diff + \
            Phtm.h_OARweight[i] * OAR_diff * OAR_diff * OAR_diff_flag;
    }

    // compare
    float mse = 0;
    for (uint i=0; i<PVCS_dose_size; i++)
    {
        float diff = h_PVCS_dose_grad[i] - h_PVCS_dose_grad_d[i];
        mse += diff * diff;
    }
    mse /= PVCS_dose_size;
    mse = sqrt(mse);

    float base = 0;
    for (uint i=0; i<PVCS_dose_size; i++)
        base += h_PVCS_dose_grad[i] * h_PVCS_dose_grad[i];
    base /= PVCS_dose_size;
    base = sqrt(base);

    cout << "The MSE value between host and device results of PVCS_dose_grad is " \
        << mse << ", while the scale of host result is " << base << endl;
    
    mse = 0;
    base = 0;
    for (uint i=0; i<PVCS_dose_size; i++)
    {
        float diff = h_element_wise_loss[i] - h_element_wise_loss_d[i];
        mse += diff * diff;
    }
    mse /= PVCS_dose_size;
    mse = sqrt(mse);

    for (uint i=0; i<PVCS_dose_size; i++)
        base += h_element_wise_loss[i] * h_element_wise_loss[i];
    base /= PVCS_dose_size;
    base = sqrt(base);
    cout << "The MSE value between host and device results of element_wise_loss is " \
        << mse << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_element_wise_loss);
    free(h_PVCS_dose_grad);
    free(h_PVCS_total_dose);
    free(h_element_wise_loss_d);
    free(h_PVCS_dose_grad_d);
    checkCudaErrors(cudaFree(d_element_wise_loss));
    checkCudaErrors(cudaFree(d_PVCS_dose_grad));
    checkCudaErrors(cudaFree(d_PVCS_total_dose));
}

void E2E::module_test_reduction()
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint size = convolved_fluence_map_dimension * convolved_fluence_map_dimension * 640;
    float* h_input = (float*)malloc(size*sizeof(float));
    float* h_loss_d = (float*)malloc(sizeof(float));
    float* d_input;
    float* d_out0;
    float* d_loss;
    checkCudaErrors(cudaMalloc((void**)&d_input, size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_out0, REDUCTION_BLOCK_SIZE*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_loss, sizeof(float)));

    uint norm = 65536;
    for (uint i=0; i<size; i++)
        h_input[i] = (float)(rand() % norm) / (norm - 1);
    
    checkCudaErrors(cudaMemcpy(d_input, h_input, size*sizeof(float), cudaMemcpyHostToDevice));

    // device compute
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "Device compute" << endl;
    cudaEventRecord(start);
    reduction(d_input, size, d_out0, d_loss, 0, 0);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << "ms" << endl;

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_loss_d, d_loss, sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    cout << "Host compute" << endl;
    double h_loss = 0;
    for (uint i=0; i<size; i++)
        h_loss += h_input[i];

    cout << "Host result: " << h_loss << ", device result: " << *h_loss_d << endl;

    // cleanup
    free(h_input);
    free(h_loss_d);
    checkCudaErrors(cudaFree(d_input));
    checkCudaErrors(cudaFree(d_out0));
    checkCudaErrors(cudaFree(d_loss));
}

void E2E::module_test_fluence_map_update(beam& Beam)
{
    uint fluence_map_size = FM_dimension * FM_dimension;
    uint extended_fluence_map_dimension = FM_dimension + 4 * FM_convolution_radius;
    uint extended_fluence_map_size = extended_fluence_map_dimension * extended_fluence_map_dimension;
    float* h_fluence_map_grad = (float*)malloc(fluence_map_size*sizeof(float));
    float* h_fluence_map = (float*)malloc(fluence_map_size*sizeof(float));
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size*sizeof(float));
    float* h_extended_fluence_map_d = (float*)malloc(extended_fluence_map_size*sizeof(float));
    float* d_squared_grad;
    float* d_final_norm;
    checkCudaErrors(cudaMalloc((void**)&d_squared_grad, fluence_map_size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_final_norm, sizeof(float)));

    // initialization
    uint norm = 65525;
    for (uint i=0; i<fluence_map_size; i++)
    {
        h_fluence_map_grad[i] = (float)(rand() % norm) / (norm - 1);
        h_fluence_map[i] = (float)(rand() % norm) / (norm - 1);
    }

    for (uint i=0; i<extended_fluence_map_size; i++)
        h_extended_fluence_map[i] = 0;
    
    for (uint i=0; i<FM_dimension; i++)
    {
        uint idx_i = i * FM_dimension;
        uint extended_idx_i = (i + 2 * FM_convolution_radius) * extended_fluence_map_dimension;
        for (uint j=0; j<FM_dimension; j++)
        {
            uint idx_j = idx_i + j;
            uint extended_idx_j = extended_idx_i + j + 2 * FM_convolution_radius;
            h_extended_fluence_map[extended_idx_j] = h_fluence_map[idx_j];
        }
    }

    // load to device
    checkCudaErrors(cudaMemcpy(Beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Beam.d_fluence_grad, h_fluence_map_grad, \
        fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    
    // device compute
    float step_size = 0.5;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cout << "Device compute" << endl;
    cudaEventRecord(start);
    Beam.fluence_map_update(0, d_final_norm, d_squared_grad, step_size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Device compute time elapse: " << milliseconds << endl;

    // offload to CPU
    checkCudaErrors(cudaMemcpy(h_extended_fluence_map_d, Beam.d_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));

    // host compute
    cout << "Host compute" << endl;
    double h_norm_final = 0;
    for (uint i=0; i<fluence_map_size; i++)
        h_norm_final += h_fluence_map_grad[i] * h_fluence_map_grad[i];
    h_norm_final /= fluence_map_size;
    h_norm_final = sqrt(h_norm_final);
    for (uint i=0; i<FM_dimension; i++)
    {
        uint idx_i = i * FM_dimension;
        uint extended_idx_i = (i + 2 * FM_convolution_radius) * extended_fluence_map_dimension;
        for (uint j=0; j<FM_dimension; j++)
        {
            uint idx_j = idx_i + j;
            uint extended_idx_j = extended_idx_i + j + 2 * FM_convolution_radius;
            h_extended_fluence_map[extended_idx_j] -= step_size * h_fluence_map_grad[idx_j] / h_norm_final;
            h_extended_fluence_map[extended_idx_j] = h_extended_fluence_map[extended_idx_j] * (h_extended_fluence_map[extended_idx_j] > 0);
        }
    }

    // compare
    float mse = 0;
    for (uint i=0; i<extended_fluence_map_size; i++)
    {
        float diff = h_extended_fluence_map[i] - h_extended_fluence_map_d[i];
        mse += diff * diff;
    }
    mse /= extended_fluence_map_size;
    mse = sqrt(mse);
    float base = 0;
    for (uint i=0; i<extended_fluence_map_size; i++)
        base += h_extended_fluence_map[i] * h_extended_fluence_map[i];
    base /= extended_fluence_map_size;

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_fluence_map_grad);
    free(h_fluence_map);
    free(h_extended_fluence_map);
    free(h_extended_fluence_map_d);
    checkCudaErrors(cudaFree(d_squared_grad));
    checkCudaErrors(cudaFree(d_final_norm));
}