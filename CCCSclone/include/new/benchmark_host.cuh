#ifndef __BENCHMARK_HOST_CUH__
#define __BENCHMARK_HOST_CUH__

#include "brain_defs.h"
#include "beam.h"
#include "kernel.h"
#include "binary_io.h"
#include <vector>

namespace dev {
    int beamsInit(std::vector<old::BEAM>& beam_arr, int verbose=0);

    float inline rand01();
    
    int radconvTextureBench(
        old::MONO_KERNELS          *mono,
        old::CONSTANTS             *constants,
        std::vector<old::BEAM>&    beams,
        int                        nrays
    );

    void radconvolveComputeBench (
        old::MONO_KERNELS *mono, old::CONSTANTS* constants, int nrays, old::BEAM& this_beam,
        old::PILLAR_GRID& hPG, int* dpg_beamletIdx, float2* dpg_beamletAngles,
        float3* dpg_pillarStartCoords, float3* dpg_beamletIsocenters,
        const std::vector<dim3>& rayGrid, const dim3& rayBlock,
        const std::vector<old::REV_DATA>& rev, float* d_fluence_map,
        const std::vector<dim3>& conGrid, const std::vector<dim3>& conBlock, 
        const std::vector<uint>& memsize, const dim3& packedGrid, const dim3& tileBlock,
        int dc
    );
}

#endif