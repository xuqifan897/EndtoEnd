#ifndef __DOSECALCIO_MACROS_H__
#define __DOSECALCIO_MACROS_H__

// CUDA device function qualifier for use in host/device compilable library code
#ifdef __CUDA_ARCH__
    #define CUDEV_FXN __host__ __device__
#else
    #define CUDEV_FXN
#endif

// Convert cuda vector types to c-style arrays
#define VECT2ARR(a, v) a[0] = v.x; a[1] = v.y;
#define VECT3ARR(a, v) a[0] = v.x; a[1] = v.y; a[2] = v.z;
// Convert c-style array to cuda vector types
#define ARR2VECT(v, a) v.x = a[0]; v.y = a[1];
#define ARR3VECT(v, a) v.x = a[0]; v.y = a[1]; v.z = a[2];

// std::cout formatting
#define FORMAT_3VEC(v) "("<<v.x<<", "<<v.y<<", "<<v.z<<")"
#define FORMAT_2VEC(v) "("<<v.x<<", "<<v.y<<")"

// use device intrinsics for device execution (faster but less accurate hardware solutions)
// TODO: -x given to intrinsic cos/sin gives NaN (investigate if true)
#if defined(__CUDA_ARCH__)
    // use device intrinsics for device execution (faster but less accurate hardware solutions)
    #define fast_cosf(x)                __cosf(x)
    #define fast_sinf(x)                __sinf(x)
    #define fast_sincosf(x, sptr, cptr) __sincosf(x, sptr, cptr)
    #define fast_powf(x, n)             __powf(x, n)
    #define fast_sq(x)                  __fmul_rn(x,x)
#else
    // fallback to helper_math.h definitions (which auto-fallback to host functions for gcc compiled code)
    #define fast_cosf(x)                cosf(x)
    #define fast_sinf(x)                sinf(x)
    #define fast_sincosf(x, sptr, cptr) sincosf(x, sptr, cptr)
    #define fast_powf(x, n)             powf(x, n)
    #define fast_sq(x)                  x*x
#endif

#endif // __DOSECALCIO_MACROS_H__