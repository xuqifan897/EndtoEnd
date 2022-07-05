#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <fstream>
#include <vector>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

extern "C"
void FCBBPVCSDoseGrad(cudaSurfaceObject_t FCBB_PVCS_dose_grad_surface, float* elementWiseLoss, \
    float* d_FCBB_PVCS_dose, float* PTV_weight, float* PTV_target, \
    float* OAR_weight, float* OAR_target, uint dimension[3], cudaStream_t stream);

void beam::calc_FCBB_PVCS_dose_grad(phantom& Phtm, float** d_elementWiseLoss, \
        float* d_PVCS_total_dose, cudaStream_t stream)
{
    /* the loss function is defined as $$PTV_weight .* ||PVCS_dose - PTV_target||_2^2 +
     OAR_weight .* ||(PVCS_dose - OAR_target)_+||_2^2$$. So the loss function should be 
     2 * PTV_weight .* (PVCS_dose - PTV_target) + OAR_weight .* (PVCS_dose - OAR_target)_+ */
    
    if (! Phtm.pitchPadding)
    {
        cout << "It is required that pitchPad() must be called" << endl;
        exit;
    }

    if (! FCBB_PVCS_dose_grad_init)
    {
        cout << "FCBBStaticInit static member function is not called!" << endl;
        exit;
    }

    if (*d_elementWiseLoss == nullptr)
    {
        cout << "d_elementWiseLoss has not been initialized" << endl;
        exit;
    }
    uint dimension[3]{Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    FCBBPVCSDoseGrad(FCBB_PVCS_dose_grad_surface, *d_elementWiseLoss, \
        d_PVCS_total_dose, Phtm.d_PTVweight, Phtm.d_PTVtarget, \
        Phtm.d_OARweight, Phtm.d_OARtarget, dimension, stream);
}

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output);

void E2E::test_calc_FCBB_PVCS_dose_grad(vector<beam>& beams, phantom& Phtm)
{
    if (! beam::FCBB_PVCS_dose_grad_init)
        beam::FCBBStaticInit(Phtm);

    beams[0].FCBBinit(Phtm);

    string inputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/doseRand.dat"};
    uint size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* rand = (float*)malloc(size*5*sizeof(float));
    ifstream inFile(inputPath);
    inFile.read((char*)rand, size*5*sizeof(float));
    inFile.close();

    float* h_FCBB_PVCS_dose = rand;
    float* h_PTVweight = rand + size;
    float* h_PTVtarget = h_PTVweight + size;
    float* h_OARweight = h_PTVtarget + size;
    float* h_OARtarget = h_OARweight + size;

    checkCudaErrors(cudaMemcpy(beams[0].d_FCBB_PVCS_dose, h_FCBB_PVCS_dose, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVweight, h_PTVweight, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVtarget, h_PTVtarget, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARweight, h_OARweight, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARtarget, h_OARtarget, size*sizeof(float), cudaMemcpyHostToDevice));

    float* d_elementWiseLoss = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_elementWiseLoss, size*sizeof(float)));
    float* h_elementWiseLoss = (float*)malloc(size*sizeof(float));

    // // for debug purposes
    // checkCudaErrors(cudaMalloc((void**)&dose_debug, size*sizeof(float)));

    beam::calc_FCBB_PVCS_dose_grad(Phtm, &d_elementWiseLoss, beams[0].d_FCBB_PVCS_dose);
    // float reduce_value = reduction(d_elementWiseLoss, size);
    // cout << reduce_value << endl;

    // // for debug purposes
    // float* h_dose_debug = (float*)malloc(size*sizeof(float));
    // checkCudaErrors(cudaMemcpy(h_dose_debug, dose_debug, size*sizeof(float), cudaMemcpyDeviceToHost));
    // string outputPath_{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_dose_grad_.dat"};
    // ofstream outFile_(outputPath_);
    // outFile_.write((char*)h_dose_debug, size*sizeof(float));
    // outFile_.close();

    float* h_FCBB_PVCS_dose_grad = (float*)malloc(size*sizeof(float));
    float* d_FCBB_PVCS_dose_grad = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_PVCS_dose_grad, size*sizeof(float)));
    
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(Phtm.dimension[0] / blockSize.x, Phtm.dimension[1] / blockSize.y, \
        Phtm.pitch / blockSize.z);
    testReadPVCSTexture(gridSize, blockSize, beam::FCBB_PVCS_dose_grad_texture, d_FCBB_PVCS_dose_grad);

    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose_grad, d_FCBB_PVCS_dose_grad, size*sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_dose_grad.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_FCBB_PVCS_dose_grad, size*sizeof(float));
    outFile.close();

    checkCudaErrors(cudaMemcpy(h_elementWiseLoss, d_elementWiseLoss, size*sizeof(float), cudaMemcpyDeviceToHost));
    outputPath = "/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_element_wise_loss.dat";
    outFile.open(outputPath);
    outFile.write((char*)h_elementWiseLoss, size*sizeof(float));
    outFile.close();
}