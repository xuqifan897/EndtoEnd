#ifndef __DEBUGLOG_H__
#define __DEBUGLOG_H__

#include "kernel.h"
#include "brain_defs.h"
#include "paths.h"

#include "boost/filesystem.hpp"

namespace old
{
    int debugLog(SHM_DATA* datavols, MONO_KERNELS* mono, CONSTANTS* constants, std::vector<BEAM>& beams);
    int datavols_log(CONSTANTS* constants, SHM_DATA* datavols);
    int mono_kernel_log(MONO_KERNELS* mono);
    int constants_log(CONSTANTS* constants);
    int beams_log(std::vector<BEAM>& beams);
    int texKern_log(CONSTANTS* constants);
    int texSpectrum_log(MONO_KERNELS* mono_kernels, const boost::filesystem::path& debugDir);
    int texDens_log(CONSTANTS* constants, const boost::filesystem::path& debugDir);
    int write_g_coords_log(const float3* d_g_coords_log, 
        CONSTANTS* constants, const boost::filesystem::path& debugDir);

    __global__ void readTexture2D(float* output, cudaTextureObject_t texObj, int width, int hight);
    __global__ void readTexture3D(float* output, cudaTextureObject_t texObj, int width, int hight, int depth);

    int logCudaVector_txt(float* d_source, int size, 
        const boost::filesystem::path& destFile);

    int logCudaVector_txt(int* d_source, int size, 
        const boost::filesystem::path& destFile);

    int logCudaVector_txt(float2* d_source, int size,
        const boost::filesystem::path& destFile);

    int logCudaVector_txt(float3* d_source, int size,
        const boost::filesystem::path& destFile);

    std::ostream& operator<<(std::ostream& os, const float3& obj);
    std::ostream& operator<<(std::ostream& os, const int3& obj);
    std::ostream& operator<<(std::ostream& os, const float2& obj);
    std::ostream& operator<<(std::ostream& os, const int2& obj);
    std::ostream& operator<<(std::ostream& os, const uint2& obj);
    std::ostream& operator<<(std::ostream& os, const dim3& obj);
}

#endif