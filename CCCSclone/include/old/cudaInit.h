#ifndef __CUDAINIT_H__
#define __CUDAINIT_H__

#include "brain_defs.h"
#include "kernel.h"
#include <cuda_runtime.h>

namespace old
{
    int initCudaConstandTex(
        SHM_DATA     *datavols,        // host-side data arrays
        MONO_KERNELS *mono,            // spectrum data
        CONSTANTS    *constants,       // calculation paramaters/information
        int          nrays,            // total number of rays
        bool         verbose=false,    // output switch
        bool         timing=false,     // output switch
        bool         debug_write=false // write debug data to file
    );

    int freeCudaTexture();

    struct DEVICE_CONV_DATA {
        /* float *terma; */
        float *dose;                    // storage of final dose for each beam
        float *device_dose;             // cleared and set for each beam calculation on each device
        float *revDens;
        float *revTerma;
        cudaArray *term_Array;
        cudaArray *dose_Array;
        cudaTextureObject_t texTerma;   // texture object created/destroyed for each beam - array mem is once initialized for all beams
        cudaTextureObject_t texDose;
        cudaSurfaceObject_t surfDose;
    };

    extern cudaArray_t d_kern_array;
    extern cudaArray_t d_spectrum_array;
    extern cudaArray_t d_dens_array;

    extern cudaTextureObject_t texSpectrum;
    extern cudaTextureObject_t texKern;
    extern cudaTextureObject_t texDens;

    extern DEVICE_CONV_DATA device_data;
    extern float* device_dose_accm;

    // variables held in each device's constant memory
    // (data init by host using cudaMemCpyToSymbol() in initCudaConstants)
    __constant__ float KERN_RADII[N_KERNEL_RADII];

    void makeTexObject(cudaTextureObject_t *texobj, cudaArray *res, int dims,
        cudaTextureAddressMode addressmode, cudaTextureFilterMode filtermode,
        bool normalizedCoords=false, 
        cudaTextureReadMode readmode=cudaTextureReadMode::cudaReadModeElementType,
        cudaResourceType restype=cudaResourceType::cudaResourceTypeLinear);

    void makeSurfObject(cudaSurfaceObject_t *surfobj, cudaArray *res);
}

#endif