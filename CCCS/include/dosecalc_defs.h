#ifndef __DOSECALC_DEFS_H__
#define __DOSECALC_DEFS_H__

//#include <cuda_runtime_api.h>  // probably don't need - NVCC should define for us
#include <helper_cuda.h>
#include <helper_math.h>
#include <math_constants.h> // CUDART_PI_F

#define PI CUDART_PI_F

// hard coded paths to data folders & files
#define RESULTS_DIR "."

// dose calculation params
#define DEFAULT_DEVICE_COUNT 4
#define MAXIMUM_DEVICE_COUNT 16  // number of concurrent GPUs to use

// BEV CCCS Convolution params
#define REV_PAD  0.0f            // expansion of REV volume by this amount on all sides [unit: cm]

// It is not beneficial to pool all threads into X-dim if TILE_DIM_X will exceed the length of a beamlet
// Do not exceed the following or occupancy will suffer:
//   CC3.0+  512 threads/block
//   CC2.0   192 threads/block
//   CC<2.0  128 threads/block
#define TILE_DIM_X 32
#define TILE_DIM_Y 4


#include "DoseCalcIO/macros.h"
#endif // __DOSECALC_DEFS_H__
