#ifndef __VECTOR_OPS_H__
#define __VECTOR_OPS_H__

#include <cuda_runtime.h>

#include "./macros.h"

CUDEV_FXN inline int3 float2int_floor(const float3& a) {
    return int3 { (int)floorf(a.x), (int)floorf(a.y), (int)floorf(a.z) };
}
CUDEV_FXN inline uint3 float2uint_floor(const float3& a) {
    return uint3 { (uint)floorf(a.x), (uint)floorf(a.y), (uint)floorf(a.z) };
}
CUDEV_FXN inline int3 float2int_ceil(const float3& a) {
    return int3 { (int)ceilf(a.x), (int)ceilf(a.y), (int)ceilf(a.z) };
}
CUDEV_FXN inline uint3 float2uint_ceil(const float3& a) {
    return uint3 { (uint)ceilf(a.x), (uint)ceilf(a.y), (uint)ceilf(a.z) };
}

// NARROWING OPS
CUDEV_FXN inline unsigned int product(const uint3& a) { return a.x * a.y * a.z; }
CUDEV_FXN inline int product(const int3& a) { return a.x * a.y * a.z; }

#endif // __VECTOR_OPS_H__
