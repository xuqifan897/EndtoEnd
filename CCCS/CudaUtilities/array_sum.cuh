#ifndef __CUDA_UTILITIES_ARRAY_SUM_CUH__
#define __CUDA_UTILITIES_ARRAY_SUM_CUH__

#include <cstdint>

// sum same-sized arrays
__global__ void cudaSum( float *sum, float *single, uint64_t size );

// accumulate smaller array into larger array with array sizes and small array start indices
__global__ void cudaSumSubArray( float *sum, float *single, uint3 smin, uint3 single_size, uint3 sum_size );

template<typename T>
__global__ void cudaResample(
    T* out,
    float3 o_start,
    uint3  o_size,
    float3 o_spacing,

    cudaTextureObject_t in,
    float3 i_start,
    uint3  i_size,
    float3 i_spacing
) {
  // Resample "in" to "out" using their respective coordinate systems
  uint3 oid = uint3{
    blockIdx.x * blockDim.x + threadIdx.x,
    blockIdx.y * blockDim.y + threadIdx.y,
    blockIdx.z * blockDim.z + threadIdx.z
  };

  float3 iid = float3{
    (oid.x*o_spacing.x + o_start.x - i_start.x)/i_spacing.x,
    (oid.y*o_spacing.y + o_start.y - i_start.y)/i_spacing.y,
    (oid.z*o_spacing.z + o_start.z - i_start.z)/i_spacing.z
  };

  if (oid.x <= o_size.x && oid.y <= o_size.y && oid.z <= o_size.z) {
    T val = tex3D<T>(in, iid.x+0.5f, iid.y+0.5f, iid.z+0.5f);
    uint olid = oid.x + o_size.x*(oid.y + o_size.y*oid.z);
    out[olid] = val;
  }
}

template<typename T>
__global__ void cudaThreshold(
    T* out,
    T* in,
    uint3    size,
    T    thresh
) {
  uint3 oid = uint3{
    blockIdx.x * blockDim.x + threadIdx.x,
    blockIdx.y * blockDim.y + threadIdx.y,
    blockIdx.z * blockDim.z + threadIdx.z
  };
  if (oid.x <= size.x && oid.y <= size.y && oid.z <= size.z) {
    uint olid = oid.x + size.x*(oid.y + size.y*oid.z);
    out[olid] = (in[olid]>=thresh)?1:0;
  }
}

#endif // __CUDA_UTILITIES_ARRAY_SUM_CUH__
