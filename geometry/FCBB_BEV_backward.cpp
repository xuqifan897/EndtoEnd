#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "geom.h"
#include "args.h"
#include "optim.h"

using namespace E2E;
using namespace std;

extern "C"
void BEVDoseBackward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    float phantom_size[3], float phantom_iso[3], \
    float* d_convolved_fluence_grad, \
    cudaTextureObject_t phantom_texture, \
    float* d_FCBB_BEV_dose_grad, \
    FCBBkernel* FCBB_kernel, \
    cudaStream_t stream);

void beam::BEV_dose_backward(phantom& Phtm, FCBBkernel* kernel, cudaStream_t stream)
{
    float phantom_size[3]{Phtm.dimension[0]*Phtm.voxelSize, Phtm.dimension[1]*\
        Phtm.voxelSize, Phtm.pitch*Phtm.voxelSize};
    float phantom_iso[3]{Phtm.isocenter[0], Phtm.isocenter[1], Phtm.isocenter[2]};
    if (this->d_convolved_fluence_map_grad == nullptr)
    {
        cout << "d_convolved_fluence_map_grad is not initialized. E2E::beams_init() not called." << endl;
        exit;
    }
    BEVDoseBackward(this->zenith, this->azimuth, this->SAD, this->pixel_size, \
        this->sampling_range[0], this->sampling_range[1], this->sampling_points, \
        phantom_size, phantom_iso, \
        this->d_convolved_fluence_map_grad, \
        Phtm.tex, \
        this->d_FCBB_BEV_dose_grad, \
        kernel, \
        stream);
}

void E2E::test_FCBB_BEV_backward(std::vector<beam>& beams, phantom& Phtm)
{
    if (!beam::FCBB_PVCS_dose_grad_init)
        beam::FCBBStaticInit(Phtm);
    
    beam& this_beam = beams[0];
    this_beam.FCBBinit(Phtm);

    // first, we conduct random number experiment
    uint convolved_fluence_map_size = this_beam.convolved_fluence_map_dimension[0] * \
        this_beam.convolved_fluence_map_dimension[1];
    float* h_FCBB_BEV_dose_grad = (float*)malloc(convolved_fluence_map_size*sizeof(float));
    srand(1008611);
    uint norm = 1024;
    for (uint i=0; i<convolved_fluence_map_size; i++)
        h_FCBB_BEV_dose_grad[i] = (float)(rand()%norm) / (norm-1);
    for (uint i=0; i<10; i++)
        cout << h_FCBB_BEV_dose_grad[i];
    cout << endl;
}