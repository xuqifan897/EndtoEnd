#include "array_sum.cuh"

// summation function across multiple streams on a single GPU
__global__ void cudaSum( float *sum, float *single, uint64_t size ) {
    std::uint64_t tid = threadIdx.x + blockDim.x*(blockIdx.x + gridDim.x*blockIdx.y);

    if ( tid >= size )
        return;

    __syncthreads();
    sum[tid] += single[tid];
}

__global__ void cudaSumSubArray ( float *sum, float *single, uint3 smin, uint3 single_size, uint3 sum_size )
{
    // tid - threadid
    // vid - linearized volume index
    unsigned int tid = threadIdx.x + blockDim.x*(blockIdx.x + gridDim.x*blockIdx.y);
    unsigned int X = ( tid % ( single_size.x * single_size.y ) ) % single_size.x;
    unsigned int Y = ( tid % ( single_size.x * single_size.y ) ) / single_size.x;
    unsigned int Z = ( tid / ( single_size.x * single_size.y ) );
    unsigned int vid = (X + smin.x) + (sum_size.x * ((Y + smin.y) + (sum_size.y * (Z + smin.z))));

    if ( tid >= (single_size.x * single_size.y * single_size.z) ||
            vid >= (sum_size.x * sum_size.y * sum_size.z) )
        return;

    __syncthreads();
    sum[vid] += single[tid];
    /* sum[tid] += single[tid]; */  // TODO Why is this here?
}
