#ifndef __CUDASIDDON_CUH__
#define __CUDASIDDON_CUH__

#include "DoseCalcIO/beam.h"
#include "server/brain_defs.h"
#include "CudaUtilities/geometry.cuh"

__device__ float cudaRaytrace_device(
        float3              this_voxel,
        float3              start,
        float3              end,
        float3              voxel,
        float3              beam_source,
        cudaTextureObject_t tex3Dvol
);

// kernel combines the radiological effective depth calculation (based on Siddon's algorithm)
// with the TERMA calculation
__device__ float cudaSiddon_device(
        float3              this_voxel,
        float*              fluence_map,
        float3              start,
        float3              end,
        float3              voxel,
        float               beamhard_correct,
        float3              beam_source,
        float3              beam_direction,
        float3              beam_isocenter,
        float               beam_sad,
        float2              beamlet_size,
        uint2               fmap_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        int                 dkerns,
        cudaTextureObject_t tex3Ddens,
        cudaTextureObject_t specTex,
        int                 beamletnum=-1
        );

__global__ void cudaSiddon(
        float*              terma,
        float*              fluence_map,
        float3              start,
        float3              f_size,
        float3              voxel,
        uint3               calc_bbox_start,
        uint3               calc_bbox_size,
        float               beamhard_correct,
        float3              beam_source,
        float3              beam_direction,
        float3              beam_isocenter,
        float               beam_sad,
        float2              beamlet_size,
        uint2               fmap_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        int                 dkerns,
        cudaTextureObject_t tex3Ddens,
        cudaTextureObject_t specTex,
        int                 beamletnum=-1
        );

#endif //__CUDASIDDON_CUH__


