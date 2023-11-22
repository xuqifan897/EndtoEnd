#ifndef __CUDADOSECALC_H__
#define __CUDADOSECALC_H__

#include "brain_defs.h"
#include "beam.h"
#include "kernel.h"
#include "binary_io.h"

#include <vector>

namespace old
{
    int radconvolveTexture (
        MONO_KERNELS        *mono,
        CONSTANTS           *constants,
        std::vector<BEAM>&  beams,
        int                 nrays,
        RES_LOG& result
    );

    int radconvolvePrep(CONSTANTS* constants, PILLAR_GRID* hPG, 
        int nrays, std::vector<REV_DATA>& rev, BEAM& this_beam, std::vector<dim3>& rayGrid, dim3& rayBlock,
        std::vector<dim3>& conGrid, std::vector<dim3>& conBlock, std::vector<uint>& memsize,
        uint3& max_actual_rev_size);
    
    int radconvolveCompute(
        MONO_KERNELS *mono, CONSTANTS* constants, int nrays, BEAM& this_beam,
        PILLAR_GRID& hPG, int* dpg_beamletIdx, float2* dpg_beamletAngles,
        float3* dpg_pillarStartCoords, float3* dpg_beamletIsocenters,
        const std::vector<dim3>& rayGrid, const dim3& rayBlock,
        const std::vector<REV_DATA>& rev, float* d_fluence_map,
        const std::vector<dim3>& conGrid, const std::vector<dim3>& conBlock, 
        const std::vector<uint>& memsize, const dim3& packedGrid, const dim3& tileBlock,
        int dc, BEAM_LOG& result);

    void update_extents(float3& currMin, float3& currMax, const float3& thisPt);

    class ArraySizeError : public virtual std::runtime_error
    {
    public:
        float theta;
        float phi;
        uint3 oldsize;
        REV_DATA* rev;
        ArraySizeError(float theta, float phi, uint3 oldsize, REV_DATA* rev) :
            std::runtime_error{""}, theta{theta}, phi{phi}, oldsize{oldsize}, rev{rev} {}
    };

    void findREV(
        REV_DATA  *rev,
        CONSTANTS *constants,
        BEAM      *beam,
        float3    lim_min,
        float3    lim_max,
        float     theta,         // azimuthal angle of the convolution ray [0->pi]
        float     phi,           // zenithal angle of the convolution ray [0->2pi]
        bool      verbose=false);
    
    __global__ void
    cudaBeamletRaytrace(
        float               *packbDens,                 // storage for BEV resampled density
        float               *packbTerm,                 // storage for BEV resampled terma
        float3              beam_source,
        float2              beamlet_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        float3              f_pg_gridDims,
        float3              f_pg_pillarDims,
        int                 pg_numBeamlets,
        int                 pg_wallThickness,
        int2                pg_numPillars,
        int*                pg_beamletIdx,
        float2*             pg_beamletAngles,
        float3*             pg_pillarStartCoords,
        float3*             pg_beamletIsocenters,
        float               kern_theta,
        float               kern_phi,                   // kernel rotation angles
        float3              revStart,                   // REV volume limit coords in XYZ coord system
        uint3               rev_size,                   // size of resampled REV data volumes
        uint3               max_rev_size,
        float3              start,
        float3              voxel,
        float3              rev_voxelsize,
        float3              f_calc_bbox_start,
        cudaTextureObject_t tex3Ddens,
        // cudaTextureObject_t tex3Dter,
        float3*             g_coords_log,

        float* d_fluence_map,
        float3 f_size,
        uint3 calc_bbox_size,
        float beamhard_correct,
        float3 beam_direction,
        float3 beam_isocenter,
        float beam_sad,
        uint2 fmap_size,
        uint dkerns,
        cudaTextureObject_t specTex,
        int requested_beamletnum = -1                   // if <0, raytrace for all beamlets
    );

    __global__ void
    PackRowConvolve(float *bDens,
        float *bTerma,
        cudaSurfaceObject_t surfDoseObj,
        float f_kern_wt_idx,
        uint3 rev_size,                     // size of resampled REV data volumes
        uint3 max_rev_size,
        float delr,
        uint nradii,
        uint ntheta,
        uint nphi,
        cudaTextureObject_t kernTex,
        float* debugProbe
    );

    __global__ void
    PackedREVtoBEVdose(
        float*              bev_packed_dose,            // beamlet-packed dose array in BEV orientation
        cudaTextureObject_t texPackedREVDose,           // packed dose array embedded in REV bounding box
        float               kern_theta, float kern_phi, // convolution direction
        float3              revStart,                    // REV volume limit coords in XYZ coord system
        float3              rev_voxelsize,
        int3                pg_gridDims
    );

    __global__ void
    UnpackBEVDosePillar (
        float*              unpacked_dose,    // output array with size matching calc_bbox for unpacked coeff. from one beamlet
        cudaTextureObject_t texPackedBEVDose, // packed dose array oriented in BEV
        float               beam_sad,
        float3              beam_source,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        float3              start,
        float3              voxel,
        float3              f_calc_bbox_start,
        uint3               calc_bbox_size,
        float3              rev_voxelsize,
        int                 pillarIndex,     // indicates which pillar to unpack
        int                 pillarIndex_x,     // indicates which pillar to unpack
        int                 pillarIndex_y,     // indicates which pillar to unpack
        int3                pg_pillarDims,
        int                 pg_wallThickness,
        int                 pg_pillarBuffer,
        float3              pillarStartCoords, // RCS coords indicating the position of the first pillar voxel
        float2              beamletAngles // beamlet divergence angles attached to this pillar
    );

    int radconvolveCompute(
        MONO_KERNELS *mono, CONSTANTS* constants, int nrays, BEAM& this_beam,
        PILLAR_GRID& hPG, int* dpg_beamletIdx, float2* dpg_beamletAngles,
        float3* dpg_pillarStartCoords, float3* dpg_beamletIsocenters,
        const std::vector<dim3>& rayGrid, const dim3& rayBlock,
        const std::vector<REV_DATA>& rev, float* d_fluence_map,
        const std::vector<dim3>& conGrid, const std::vector<dim3>& conBlock, 
        const std::vector<uint>& memsize, const dim3& packedGrid, const dim3& tileBlock,
        bool unpack2Patient, RES_LOG& beamResult);
}

#endif