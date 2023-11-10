#ifndef __DEBUGLOG_H__
#define __DEBUGLOG_H__

#include "kernel.h"
#include "brain_defs.h"
#include "paths.h"

namespace old
{
    int debugLog(SHM_DATA* datavols, MONO_KERNELS* mono, CONSTANTS* constants, std::vector<BEAM>& beams);
    int datavols_log(CONSTANTS* constants, SHM_DATA* datavols);
    int mono_kernel_log(MONO_KERNELS* mono);
    int constants_log(CONSTANTS* constants);
    int beams_log(std::vector<BEAM>& beams);
    int texKern_log(CONSTANTS* constants);
    int texSpectrum_log(MONO_KERNELS* mono_kernels);
    int texDens_log(CONSTANTS* constants);
    __global__ void readTexture2D(float* output, cudaTextureObject_t texObj, int width, int hight);
    __global__ void readTexture3D(float* output, cudaTextureObject_t texObj, int width, int hight, int depth);
}

#endif