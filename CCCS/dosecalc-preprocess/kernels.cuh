#include "dosecalc_defs.h"
#include "DoseCalcAlgo/cudaSiddon.cuh"
#include "../CudaUtilities/geometry.cuh"
#include <assert.h>

#define MAX_TEXTURE_DIM 384.f

float *iso_matrix;
cudaArray *imgArray=0;
texture<float, 3, cudaReadModeElementType> texImg;

// each thread controls a voxel and applies a simple conversion between HU and density
__global__ void
deviceHU2dens( float *CTunits, float *dens, uint3 count )
{
    int pos = threadIdx.x + blockDim.x*blockIdx.x;
    if ( pos >= count.x*count.y*count.z)
        return;

    float input = CTunits[pos];
    float temp = ( (float)input + 1024.0f ) / 1024.0f;

    if (temp < 0.05f)
        temp = 0.05f;

    dens[pos] = temp;
}

// the data volume is held in a texture
// this function resamples the data at regular intervals using texture memory's intrinsic linear interpolation
// and outputs an isotropic volume
__global__ void
cudaMakeIsotropicWithLUT( float *iso, float3 voxelSize, float iso_voxel, uint3 iso_size, float* lut_hunits, float* lut_massdens, int nlut)
{
    extern __shared__ float s[];
    // find the overall ID number of this thread, tid (thread index)
    unsigned int bidx = blockIdx.x + gridDim.x * ( blockIdx.y + gridDim.y * blockIdx.z );
    unsigned int tidx = threadIdx.x + blockDim.x * bidx;
    float* s_lut_hunits = s;
    float* s_lut_massdens = &s[nlut];

    if (blockDim.x >= nlut) {
        if (threadIdx.x < nlut) {
            s_lut_hunits[threadIdx.x] = lut_hunits[threadIdx.x];
            s_lut_massdens[threadIdx.x] = lut_massdens[threadIdx.x];
        }
    } else {
        for (int jj=0; jj<ceilf((float)nlut/blockDim.x); jj++) {
            int ii = jj*blockDim.x + threadIdx.x;
            if (ii < nlut) {
                s_lut_hunits[ii] = lut_hunits[ii];
                s_lut_massdens[ii] = lut_massdens[ii];
            }
        }
    }
    __syncthreads();

	// convert tid into 3D coordinates based on the size of the isotropic volume
    unsigned int slicesize = (iso_size.x*iso_size.y);
    unsigned int X = (tidx % slicesize) % iso_size.x;
    unsigned int Y = (tidx % slicesize) / iso_size.x;
    unsigned int Z = (tidx / slicesize);
    if (X >= iso_size.x || Y >= iso_size.y || Z >= iso_size.z) { return; }

	// convert that to physical distance in mm
    float3 pos = make_float3(X, Y, Z) * make_float3(iso_voxel);

    // convert from physical point within the isotropic volume
    // to a set of coordinates in the original data
    pos /= voxelSize;

    // same the point from the texture memory
    // 0.5f is added to each coordinate as a shift to the center of the voxel
    float value = tex3D( texImg, pos.x + 0.5f, pos.y + 0.5f, pos.z + 0.5f);

    // find low
    if (value <= s_lut_hunits[0]) {
        // lower bound
        value = s_lut_massdens[0];
    }
    else if (value >= s_lut_hunits[nlut-1]) {
        // extrapolate (just based on fit from last 2 data points)
        value = s_lut_massdens[nlut-2] +
            (value - s_lut_hunits[nlut-2]) *
            (s_lut_massdens[nlut-1]-s_lut_massdens[nlut-2])/(s_lut_hunits[nlut-1]-s_lut_hunits[nlut-2]);
    }
    else {
        // interpolate
        int lowidx = 0;
        while(s_lut_hunits[lowidx]<value && lowidx < nlut-1) {lowidx++;}
        lowidx--;
        value = s_lut_massdens[lowidx] +
            (value - s_lut_hunits[lowidx]) *
            (s_lut_massdens[lowidx+1]-s_lut_massdens[lowidx])/(s_lut_hunits[lowidx+1]-s_lut_hunits[lowidx]);
    }

	// write to output in DCS
    iso[tidx] = value;
}
__global__ void
cudaMakeIsotropic( float *iso, float3 voxelSize, float iso_voxel, uint3 iso_size )
{
    // find the overall ID number of this thread, tid (thread index)
    unsigned int bidx = blockIdx.x + gridDim.x * ( blockIdx.y + gridDim.y * blockIdx.z );
    unsigned int tidx = threadIdx.x + blockDim.x * bidx;

	// convert tid into 3D coordinates based on the size of the isotropic volume
    unsigned int slicesize = (iso_size.x*iso_size.y);
    unsigned int X = (tidx % slicesize) % iso_size.x;
    unsigned int Y = (tidx % slicesize) / iso_size.x;
    unsigned int Z = (tidx / slicesize);
    if (X >= iso_size.x || Y >= iso_size.y || Z >= iso_size.z) { return; }

	// convert that to physical distance in mm
    float3 pos = make_float3(X, Y, Z) * make_float3(iso_voxel);

    // convert from physical point within the isotropic volume
    // to a set of coordinates in the original data
    pos /= voxelSize;

    // same the point from the texture memory
    // 0.5f is added to each coordinate as a shift to the center of the voxel
    float value = tex3D( texImg, pos.x + 0.5f, pos.y + 0.5f, pos.z + 0.5f);

	// write to output in DCS
    iso[tidx] = value;
}

// This function is modified from Siddon's algorithm for a 2D accumulated projection
// note: everything is in units of centimeters
#define subrayidx   (threadIdx.x)
#define subraycount (blockDim.x)
#define s_index(ii) threadIdx.y+blockDim.y*(threadIdx.z+blockDim.z*ii)
__global__ void
cudaRaycast(
        float* fluence,
        float3 vol_start,
        uint3  vol_size,
        float3 voxel_size,
        float3 source,
        float3 isocenter,
        uint2  fmap_size,
        float2 beamlet_size,
        float  azimuth,
        float  zenith,
        float  coll,
        cudaTextureObject_t texRay
) {
	// dynamically allocated shared memory
    extern __shared__ char s_intersect[];
    s_intersect[s_index(subrayidx)] = 0;
    __syncthreads();

    int ray_X = threadIdx.y + blockIdx.y*blockDim.y;
    int ray_Z = threadIdx.z + blockIdx.z*blockDim.z;

	// threads are launched for each pixel in the 2D output map
    // ray_X/ray_Z in Fluence Coord Sys (FCS) which matches RCS before rotation
    if (ray_X >= fmap_size.x || ray_Z >= fmap_size.y) { return; }

    // coords of start and end in RCS
    float3 vol_end;
    vol_end.x = vol_start.x + __int2float_rn(vol_size.x) * voxel_size.x;
    vol_end.y = vol_start.y + __int2float_rn(vol_size.y) * voxel_size.y;
    vol_end.z = vol_start.z + __int2float_rn(vol_size.z) * voxel_size.z;

    // shift for sub-ray [-1, 0, +1]
    int subray_X = 0;
    int subray_Z = 0;
    if (subraycount > 1) {
        subray_X = (subrayidx % 3) - 1;
        subray_Z = (subrayidx / 3) - 1;
    }

	// the end point of each ray is found from the 2D coordinates on the fluence map
    // center coord of fluence map is isocenter
    // bixel coords defined in ray_X-ray_Z plane (FCS) then rotated with beam angles into RCS
    // dont let x/y in int2 scare you, it is really x/z
    // we directly define bixel coords in DCS to be consistent with storage of texRay
    float3 bixel_ctr_FCS = make_float3(
            beamlet_size.x*(ray_X + 0.5f*subray_X - fmap_size.x*0.5f + 0.5f),
            0,
            beamlet_size.y*(ray_Z + 0.5f*subray_Z - fmap_size.y*0.5f + 0.5f)
            );

    float3 bixel_ctr = inverseRotateBeamAtOriginRHS(bixel_ctr_FCS, azimuth, zenith, coll);
    bixel_ctr += isocenter;
	// avoid divisions by 0
    bixel_ctr += 1e-9f;

    // extend end of raytrace beyond fluence map plane
    float3 shortdiff = bixel_ctr - source;
    float3 sink = source + 3.f*shortdiff;

	// the vector of projection for this thread, extended completely through volume
    float3 diff = sink - source;

    float rpl = cudaRaytrace_device(
            sink,
            vol_start,
            vol_end,
            voxel_size,
            source,
            texRay
            );

	// any ray that passes through the contour area will have accumulated a value greater than 0

    if (rpl > 0.f) { s_intersect[s_index(subrayidx)] = 1; }
    __syncthreads();
    if (subrayidx == 0)  {
        float f = 0.f;
        for (int ii=0; ii<subraycount; ii++) {
            f += s_intersect[s_index(ii)];
        }

        // write out the fluence map
        fluence[ray_X + fmap_size.x * ray_Z] = float(f>0);
    }
}
