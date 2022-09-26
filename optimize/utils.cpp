#include <vector>
#include <math.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

extern "C"
void doseSum(float* result, float** sources, uint num_beams, uint size, cudaStream_t stream);

void E2E::dose_sum(std::vector<beam>& beams, phantom& Phtm, float** d_result, \
    float** d_sources, cudaStream_t stream)
{   
    uint size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    if (*d_result == nullptr)
    {
        cout << "d_result has not been initialized" << endl;
        exit;
    }
    doseSum(*d_result, d_sources, beams.size(), size, stream);
}

extern "C"
void Reduction(dim3 gridSize, dim3 blockSize, float* out, \
    float* source, uint size, uint idx, cudaStream_t stream);

void E2E::reduction(float* d_source, uint size, float* d_out0, \
    float* loss, uint idx, cudaStream_t stream)
{
    dim3 blockSize(REDUCTION_BLOCK_SIZE);
    dim3 gridSize0(REDUCTION_BLOCK_SIZE);
    dim3 gridSize1(1);

    Reduction(gridSize0, blockSize, d_out0, d_source, size, 0, stream);
    Reduction(gridSize1, blockSize, loss, d_out0, REDUCTION_BLOCK_SIZE, idx, stream);
}


void E2E::reduction_small(float* d_source, uint size, float* loss, \
    uint idx, cudaStream_t stream)
{
    dim3 blockSize(REDUCTION_BLOCK_SIZE);
    dim3 gridSize(1);
    Reduction(gridSize, blockSize, loss, d_source, size, idx, stream);
}


void E2E::test_dose_sum(std::vector<beam>& beams, phantom& Phtm)
{
    if (beams.size() != 5)
    {
        cout << "Here we assume the number of beams is 5" << endl;
        exit;
    }

    uint size_per_beam = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    string path{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/doseRand.dat"};
    float* h_doseRand = (float*)malloc(size_per_beam * beams.size() * sizeof(float));
    ifstream inFile(path);
    inFile.read((char*)h_doseRand, size_per_beam * beams.size() * sizeof(float));
    inFile.close();
    
    for (int i=0; i<beams.size(); i++)
    {
        if (beams[i].d_FCBB_PVCS_dose == nullptr)
            checkCudaErrors(cudaMalloc((void**)(&(beams[i].d_FCBB_PVCS_dose)), size_per_beam*sizeof(float)));
        checkCudaErrors(cudaMemcpy(beams[i].d_FCBB_PVCS_dose, h_doseRand + i * size_per_beam, \
            size_per_beam * sizeof(float), cudaMemcpyHostToDevice));
    }
    
    float* d_result = nullptr;
    float** h_sources = (float**)malloc(beams.size()*sizeof(float*));
    float** d_sources = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_result, size_per_beam*sizeof(float)));
    for (int i=0; i<beams.size(); i++)
        h_sources[i] = beams[i].d_FCBB_PVCS_dose;
    checkCudaErrors(cudaMalloc((void***)&d_sources, beams.size()*sizeof(float*)));
    checkCudaErrors(cudaMemcpy(d_sources, h_sources, beams.size()*sizeof(float*), cudaMemcpyHostToDevice));
    
    dose_sum(beams, Phtm, &d_result, d_sources);

    float* h_result = (float*)malloc(size_per_beam * sizeof(float));
    checkCudaErrors(cudaMemcpy(h_result, d_result, size_per_beam * sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/doseRandResult.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_result, size_per_beam*sizeof(float));
    outFile.close();
}

extern "C"
void elementWiseSquare(dim3 gridSize, dim3 blockSize, float* d_output, \
    float* d_input, uint size, cudaStream_t stream);

void E2E::element_wise_square(float* d_output, float* d_input, uint size, cudaStream_t stream)
{
    dim3 blockSize(256);
    dim3 gridSize((uint)ceil((float)size / blockSize.x));
    elementWiseSquare(gridSize, blockSize, d_output, d_input, size, stream);
}

void E2E::test_element_wise_square()
{
    uint size = 365;
    uint mod = 1024;
    float* h_array = (float*)malloc(size*sizeof(float));
    float* d_array;
    float* h_output = (float*)malloc(size*sizeof(float));
    float* h_output_host = (float*)malloc(size*sizeof(float));
    float* d_output;
    checkCudaErrors(cudaMalloc((void**)&d_array, size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_output, size*sizeof(float)));
    srand(1008611);
    for (uint i=0; i<size; i++)
    {
        h_array[i] = (float)(rand()%mod) / mod;
        h_output_host[i] = h_array[i] * h_array[i];
    }
    checkCudaErrors(cudaMemcpy(d_array, h_array, size*sizeof(float), cudaMemcpyHostToDevice));

    element_wise_square(d_output, d_array, size);

    checkCudaErrors(cudaMemcpy(h_output, d_output, size*sizeof(float), cudaMemcpyDeviceToHost));

    float diff = 0;
    for (uint i=0; i<size; i++)
        diff += abs(h_output[i] - h_output_host[i]);
    cout << diff << endl;
}

extern "C"
void fluenceMapUpdate(uint fluence_map_dimension, uint extended_fluence_map_dimension, uint extended_fluence_map_offset, \
    uint idx, float* d_extended_fluence_map, float* d_fluence_map_grad, float* d_norm, float step_size, \
    cudaStream_t stream);

void beam::fluence_map_update(uint idx, float* d_norm_final, \
    float* d_squared_grad, float step_size, cudaStream_t stream)
{
    // idx is the index of the beam, d_norm is device array containing the norm
    // d_squared_grad is the device array containing element_size squared_grad
    // d_norm is the device array containing the summed square
    if (d_squared_grad == nullptr)
    {
        cout << "d_squared_grad is not initialized!" << endl;
        exit;
    }
    uint fluence_map_size = this->fluence_map_dimension[0] * this->fluence_map_dimension[1];
    element_wise_square(d_squared_grad, this->d_fluence_grad, fluence_map_size, stream);
    
    dim3 blockSize(REDUCTION_BLOCK_SIZE);
    dim3 gridSize(1);
    Reduction(gridSize, blockSize, d_norm_final, d_squared_grad, fluence_map_size, idx, stream);
    fluenceMapUpdate(this->fluence_map_dimension[0], this->extended_fluence_map_dimension[0], 2*FM_convolution_radius, \
        idx, this->d_extended_fluence_map, this->d_fluence_grad, d_norm_final, step_size, stream);
}

void E2E::test_fluence_map_update(vector<beam>& beams)
{
    uint idx = 2;
    float step_size = 0.314;
    beam& this_beam = beams[idx];
    uint fluence_map_size = this_beam.fluence_map_dimension[0] * this_beam.fluence_map_dimension[1];
    uint extended_fluence_map_size = this_beam.extended_fluence_map_dimension[0] * this_beam.extended_fluence_map_dimension[1];

    float* d_norm_final;
    float* d_squared_grad;
    checkCudaErrors(cudaMalloc((void**)&d_norm_final, beams.size()*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_squared_grad, fluence_map_size*sizeof(float)));

    // float* h_norm_final = (float*)malloc(beams.size()*sizeof(float));
    float* h_fluence_map_grad = (float*)malloc(fluence_map_size*sizeof(float));
    float* h_fluence_map = (float*)malloc(fluence_map_size*sizeof(float));
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size*sizeof(float));

    // random initialization
    srand(1008611);
    uint norm = 4096;
    for (uint i=0; i<fluence_map_size; i++)
    {
        h_fluence_map[i] = (float)(rand() % norm) / (norm-1);
        h_fluence_map_grad[i] = (float)(rand() % norm) / (norm-1);
    }
    for (uint i=0; i<extended_fluence_map_size; i++)
        h_extended_fluence_map[i] = 0;
    for (uint i=0; i<this_beam.fluence_map_dimension[0]; i++)
    {
        for (uint j=0; j<this_beam.fluence_map_dimension[1]; j++)
        {
            uint fluence_map_idx = i * this_beam.fluence_map_dimension[1] + j;
            uint extended_fluence_map_idx = (i + 2 * FM_convolution_radius) * \
                this_beam.extended_fluence_map_dimension[1] + j + 2 * FM_convolution_radius;
            h_extended_fluence_map[extended_fluence_map_idx] = h_fluence_map[fluence_map_idx];
        }
    }

    // load to device
    checkCudaErrors(cudaMemcpy(this_beam.d_fluence_grad, h_fluence_map_grad, \
        fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(this_beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    
    // device execute
    this_beam.fluence_map_update(idx, d_norm_final, d_squared_grad, step_size);

    // offload result
    float* h_updated_extended_fluence_map = (float*)malloc(extended_fluence_map_size*sizeof(float));
    checkCudaErrors(cudaMemcpy(h_updated_extended_fluence_map, this_beam.d_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
    
    // calculate the result on CPU
    float* host_result = (float*)malloc(extended_fluence_map_size*sizeof(float));
    float Norm = 0;
    for (uint i=0; i<fluence_map_size; i++)
        Norm += h_fluence_map_grad[i] * h_fluence_map_grad[i];
    Norm /= fluence_map_size;
    Norm = sqrt(Norm);
    for (uint i=0; i<extended_fluence_map_size; i++)
        host_result[i] = 0;
    for (uint i=0; i<this_beam.fluence_map_dimension[0]; i++)
    {
        for (uint j=0; j<this_beam.fluence_map_dimension[1]; j++)
        {
            uint fluence_map_idx = i * this_beam.fluence_map_dimension[1] + j;
            uint extended_fluence_map_idx = (i + 2 * FM_convolution_radius) * \
                this_beam.extended_fluence_map_dimension[1] + j + 2 * FM_convolution_radius;
            host_result[extended_fluence_map_idx] = h_extended_fluence_map[extended_fluence_map_idx] - \
                h_fluence_map_grad[fluence_map_idx] * step_size / Norm;
        }
    }

    // // compare the difference
    // float diff = 0;
    // for (uint i=0; i<extended_fluence_map_size; i++)
    //     diff += abs(host_result[i] - h_updated_extended_fluence_map[i]);
    // diff /= extended_fluence_map_size;
    // cout << "The difference between host code and device code is " << diff << endl;

    string output_file{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/updated_host.dat"};
    ofstream outFile(output_file);
    outFile.write((char*)host_result, extended_fluence_map_size*sizeof(float));
    outFile.close();
    
    output_file = "/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/updated_device.dat";
    outFile.open(output_file);
    outFile.write((char*)h_updated_extended_fluence_map, extended_fluence_map_size*sizeof(float));
    outFile.close();

    output_file = "/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/original_host.dat";
    outFile.open(output_file);
    outFile.write((char*)h_extended_fluence_map, extended_fluence_map_size*sizeof(float));
    outFile.close();
}

extern "C"
void smoothnessCalc(float* d_extended_fluence_map, float* d_element_wise_loss, float* d_fluence_grad, float eta, cudaStream_t stream);

void beam::smoothness_calc(float eta, cudaStream_t stream)
{
    smoothnessCalc(this->d_extended_fluence_map, this->d_element_wise_fluence_smoothness_loss, \
        this->d_fluence_grad, eta, stream);
}