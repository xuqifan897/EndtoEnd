#ifndef __CUDA_UTILITIES_MACROS_H__
#define __CUDA_UTILITIES_MACROS_H__

// CUDA device function qualifier for use in host/device compilable library code
#ifdef __CUDA_ARCH__
    #define CUDEV_FXN __host__ __device__
#else
    #define CUDEV_FXN
#endif

#endif // __CUDA_UTILITIES_MACROS_H__

