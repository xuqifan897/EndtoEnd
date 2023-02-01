#ifndef __NVBBRAYCONVOLVE_DEVICE_CUH__
#define __NVBBRAYCONVOLVE_DEVICE_CUH__

#include "helper_math.h"
#include "dosecalc_defs.h"
#include "server/brain_defs.h"
#include "CudaUtilities/geometry.cuh" // coordinate rotations
#include "DoseCalcIO/kernel.h"
#include <assert.h>
#include "CudaUtilities/dev_intrinsics.cuh"

// variables held in each device's constant memory
// (data init by host using cudaMemCpyToSymbol() in initCudaConstants)
__constant__ float             KERN_RADII[N_KERNEL_RADII];

#endif // __NVBBRAYCONVOLVE_DEVICE_CUH__
