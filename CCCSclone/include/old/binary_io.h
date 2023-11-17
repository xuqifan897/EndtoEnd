#ifndef __BINARY_IO_H__
#define __BINARY_IO_H__

#include "brain_defs.h"

#include <vector>

namespace old
{
    int load_data(CONSTANTS* host, SHM_DATA* data);
    int load_density(SHM_DATA* data);
    int write_result(const std::vector<std::vector<std::vector<float>>>& result);

    struct PILLAR_GRID{
        int numBeamlets;               // number of active pillars in this fluence map
        int2 numPillars;               // number of pillars along x and z axis in pillarMap.
        int3 gridDims;                 // precalculated numPillars * numVoxelPerPillar (total voxel dims of packed array)
        int3 pillarDims;               // dimension of pillars along x, y, z directions in voxels
        float3 gridSize;               // dimensions of grid in common coordinate unit
        float3 pillarSize;             // dimensions of pillars in common coordinate unit
        float2 max_beamlet_size;       // maximum size of beamlet due to divergence/magnification
        int wallThickness = 2;         // number of REV voxels composing separating walls; 
                                       // shoule be >=2 for best results
        int pillarBuffer = 2;          // number of REV voxels buffering usable pillar space from walls 
                                       // (garbage ends up here from interpolation near walls);
                                       // should be >=2 for best results
        // arrays
        std::vector<int> beamletIdx;               // indices of actual beamlets
        std::vector<float3> pillarStartCoords;     // coords of pillar box start
        std::vector<float2> beamletAngles;         // sum of beam angle and beamlet divergence angle (.x: azimuth, .y: zenith)
        std::vector<float3> beamletIsocenters;     // RCS coords of each beamlet isocenter 
                                       // (intersect of beamlet central axis and fluence map plane)
        std::vector<float3> beamletRotationAnchors;// RCS-aligned offset from pillarStartCoords that 
                                       // defines the beamletRotation angle point
        
        // convenience functions
        __host__ __device__ inline int pillar_nvoxels() const {return pillarDims.x * pillarDims.y * pillarDims.z;}
        __host__ __device__ inline int pillar_grid_nvoxels() const {return gridDims.x * gridDims.y * gridDims.z;}

    };

    struct REV_DATA {
        float3 min_coords {}; // coords of first voxel in BEV
        float3 max_coords {}; // coords of last voxel in BEV
        uint3  size {};       // dimensions of this BEV volume
    };
}

#endif