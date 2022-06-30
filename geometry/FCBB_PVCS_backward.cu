#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>
#include "geom.h"

using namespace E2E;
using namespace std;

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