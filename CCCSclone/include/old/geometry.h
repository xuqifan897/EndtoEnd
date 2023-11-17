#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <helper_math.h>
#include "macros.h"

namespace old
{
    // Rotate vec around center using arbitrary axis
    float3 rotateAroundAxisRHS(
        const float3& p, const float3& q, const float3& r, const float& t);

    float3 rotateAroundAxisAtOriginRHS(
        const float3& p, const float3& r, const float& t);

    // convert RCS coords to BEV coords
    float3 rotateBeamRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi, const float& coll);

    float3 rotateBeamAtOriginRHS(
        const float3& vec, const float& theta, const float& phi, const float& coll);

    // convert BEV coords to RCS coords
    float3 inverseRotateBeamRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi, const float& coll);

    float3 inverseRotateBeamAtOriginRHS(
        const float3& vec, const float& theta, const float& phi, const float& coll);

    // convert bev (beamlet) coords to BEV (beam) coords
    float3 rotateBeamletRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi);

    float3 rotateBeamletAtOriginRHS(
        const float3& vec, const float& theta, const float& phi);

    // convert BEV (beam) coords to bev (beamlet) coords
    float3 inverseRotateBeamletRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi);

    float3 inverseRotateBeamletAtOriginRHS(
        const float3& vec, const float& theta, const float& phi);

    // convert BEV coords to REV coords
    float3 rotateKernelRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi);
    
    float3 rotateKernelAtOriginRHS(
        const float3& vec, const float& theta, const float& phi);

    // convert REV coords to BEV coords
    float3 inverseRotateKernelRHS(
        const float3& vec, const float3& center, const float& theta, const float& phi);

    float3 inverseRotateKernelAtOriginRHS(
        const float3& vec, const float& theta, const float& phi);
    
    
    void calcBeamletAnchors(
        float3&      start,            // out: start_anchor_coords
        float3&      end,              // out: end_anchor_coords
        float2&      beamletAngles,    // out: beam angles + beamlet divergence angles (.x: azi, .y: .zen)
        float3&      beamletIso,       // out: beamlet isocenter coords
        float3       src,              // beam src coords
        float3       iso,              // beam iso coords
        unsigned int beamlet_idx,      // fluence map linearized index
        float2       beamlet_size,     // [unit: cm]
        uint2        fmap_dims,        // number of beamlets along each axis
        float3       voxelsize,        // [unit: cm]
        float3       density_start,    // coords of start of density volume
        uint3        calc_bbox_start,  // nvoxel offset from density start
        uint3        calc_bbox_size,   // nvoxel size from calc_bbox_start
        float        gantry_angle_rad, // [unit: rad]
        float        couch_angle_rad,  // [unit: rad]
        float        coll_angle_rad,   // [unit: rad]
        int          verbose=false);
}

#endif