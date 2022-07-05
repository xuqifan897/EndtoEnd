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