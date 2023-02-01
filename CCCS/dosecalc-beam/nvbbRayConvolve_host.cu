#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstring>

// cuda multithreading helpers
#include "Utilities/multithreading.h"

#include "nvbbRayConvolve_device.cu" // device memory pointers & kernel functions
#include "dosecalc_defs.h"  // MAXIMUM_DEVICE_COUNT, MAXIMUM_STREAM_COUNT
#include "server/brain_defs.h" // CONSTANTS, SHM_DATA ...

#include "CudaUtilities/cuda_timing.cuh"
#include "CudaUtilities/manage_gpu.cuh"
#include "CudaUtilities/make_tex_surf.cuh"
#include "CudaUtilities/geometry.cuh" // coordinate rotations
#include "CudaUtilities/array_sum.cuh"
#include "DoseCalcIO/dosecalcio.h"
#include "DoseCalcAlgo/cudaSiddon.cuh"


extern CONSTANTS*        constants;
extern bool extra_verbose;
extern bool extra_debug;

// cuda Arrays for texture mapping (pointers to device global mem)
// Constant for all beams; allocated on each device for efficient access (trade speed for memory)
cudaArray* d_kern_array     [ MAXIMUM_DEVICE_COUNT ] = {nullptr};
cudaArray* d_spectrum_array [ MAXIMUM_DEVICE_COUNT ] = {nullptr};
cudaArray* d_dens_array     [ MAXIMUM_DEVICE_COUNT ] = {nullptr};

// texture reference definitions (only accessed from this module, static)
// Need one ref for each device, referencing the copy of array data on that device
static cudaTextureObject_t texSpectrum [ MAXIMUM_DEVICE_COUNT ] = {0};
static cudaTextureObject_t texKern     [ MAXIMUM_DEVICE_COUNT ] = {0};
static cudaTextureObject_t texDens     [ MAXIMUM_DEVICE_COUNT ] = {0};

// all GPU arrays needed for dose convolution - stored per device
struct DEVICE_CONV_DATA {
    float *terma;
    float *dose;                    // storage of final dose for each beam (handled in sequence)
    float *device_dose;             // cleared and set for each beam calculation on each device
    float *bevDens;
    float *bevTerma;
    cudaArray *term_Array;
    cudaArray *dose_Array;
    cudaTextureObject_t texTerma;   // texture object created/destroyed for each beam - array mem is once initialized for all beams
    cudaTextureObject_t texDose;    //
    cudaSurfaceObject_t surfDose;
};

// Unique for each device
DEVICE_CONV_DATA  device_data       [ MAXIMUM_DEVICE_COUNT ] = {0};       // holds cuda execution device specific data

// finding the extents for active calculation volume in beam's eye view (BEV) for each convolution angle
// by taking the 8 corners of the active calculation volume and rotating them into a coordinate system aligned along the
// convolution direction. New resolution in BEV defined (DEFS.H) by the following variables:
/* Get RCS coordinates of start/end voxels in REV coordinate system rotated to REV orientation
 * This is critical for converting between REV indices to RCS coords then RCS indices to sample terma/density in RCS
 */
static inline void update_extents(float3& currMin, float3& currMax, const float3& thisPt) {
    if (thisPt.x < currMin.x) currMin.x = thisPt.x;
    if (thisPt.y < currMin.y) currMin.y = thisPt.y;
    if (thisPt.z < currMin.z) currMin.z = thisPt.z;
    if (thisPt.x > currMax.x) currMax.x = thisPt.x;
    if (thisPt.y > currMax.y) currMax.y = thisPt.y;
    if (thisPt.z > currMax.z) currMax.z = thisPt.z;
    return;
}
class ArraySizeError : public virtual std::runtime_error {
    public:
        ArraySizeError(const std::string& str) : std::runtime_error(str) {}
};
static void findREV(
        REV_DATA  *rev,
        CONSTANTS *constants,
        BEAM      *beam,
        float3    lim_min,
        float3    lim_max,
        float     theta,         // azimuthal angle of the convolution ray [0->pi]
        float     phi,           // zenithal angle of the convolution ray [0->2pi]
        bool      verbose=false
) {
    // diff in RCS coords relative to origin at min
    // rotated into BEV space relative to origin at min - used to get corner coords in BEV space of terma bounding box
    /* float3 b_diff = rotateBeamAtOriginRHS(lim_max-lim_min, beam->azimuth, beam->zenith, beam->coll); */
    float3 b_diff = lim_max-lim_min;
    float3 rev_min = make_float3(std::numeric_limits<float>::max());
    float3 rev_max = make_float3(std::numeric_limits<float>::min());

    // for 8 corners of the data volume defined by RCS coordinates of BEV box, find position in REV coord sys
    for (int xoff=0; xoff<2; xoff++) {
        for (int yoff=0; yoff<2; yoff++) {
            for (int zoff=0; zoff<2; zoff++) {
                float3 input = make_float3(b_diff.x*xoff, b_diff.y*yoff, b_diff.z*zoff);
                // apply rotation for the convolution ray and evaluate extents in REV space
                float3 output = rotateKernelAtOriginRHS(input, theta, phi);
                // set beam's eye view extents
                update_extents(rev_min, rev_max, output);
    } } }

    // PAD REV volume
    rev_min -= REV_PAD;
    rev_max += REV_PAD;


    // NVB calculation grid spacing is set by constants->rev_{lat,long}spacing so that conv. step size (longspacing)
    // can be independently adjusted
    // remember: Pillar_grid to REV orientation is XYZ -> ZXY and rev size reflects REV orientation
    rev->size = make_uint3(
            static_cast<unsigned int>( ceil(fabsf(rev_max.y-rev_min.y) / constants->rev_longspacing) ),
            static_cast<unsigned int>( ceil(fabsf(rev_max.z-rev_min.z) / constants->rev_latspacing) ),
            static_cast<unsigned int>( ceil(fabsf(rev_max.x-rev_min.x) / constants->rev_latspacing) ) );

    // shrink this bev
    uint3 oldsize = rev->size;
    float3 adjust = {};
    bool was_shrunk = false;
    if (rev->size.x > constants->max_rev_size.x) {
        adjust.x = constants->rev_longspacing*((int)rev->size.x - (int)constants->max_rev_size.x);
        rev->size.x = constants->max_rev_size.x;
        was_shrunk = true;
    }
    if (rev->size.y > constants->max_rev_size.y) {
        adjust.y = constants->rev_latspacing*((int)rev->size.y - (int)constants->max_rev_size.y);
        rev->size.y = constants->max_rev_size.y;
        was_shrunk = true;
    }
    if (rev->size.z > constants->max_rev_size.z) {
        adjust.z = constants->rev_latspacing*((int)rev->size.z - (int)constants->max_rev_size.z);
        rev->size.z = constants->max_rev_size.z;
        was_shrunk = true;
    }
    if (was_shrunk) { rev_max -= make_float3(adjust.z, adjust.x, adjust.y); }

    // store limits in REV coordinates relative to rotated lim_min, lim_max
    rev->min_coords = rev_min;
    rev->max_coords = rev_max;

    if (verbose) {
        printf(" Theta               :  %5.1f deg\n"            , theta*180/PI );
        printf(" Phi                 :  %5.1f deg\n"            , phi*180/PI );
        printf(" REV->size (PG:YZX)  :    %5d x   %5d x   %5d\n", rev->size.x, rev->size.y, rev->size.z);
        printf(" pbev_min            :  %7.2f x %7.2f x %7.2f\n", rev->min_coords.x, rev->min_coords.y, rev->min_coords.z);
        printf(" pbev_max            :  %7.2f x %7.2f x %7.2f\n", rev->max_coords.x, rev->max_coords.y, rev->max_coords.z);
    }
    if (was_shrunk) {
        std::ostringstream msg;
        msg << set_color(COLOR::RED)<<"Convolution array size for ray (az:"<<theta*180.f/PI<<"|zen:"<<phi*180.f/PI<<" deg) was reduced"
            << " from ("<<oldsize.x<<","<<oldsize.y<<","<<oldsize.z<<") to ("<<rev->size.x<<","<<rev->size.y<<","<<rev->size.z<<"). "
            << "Results are likely to contain errors"<< set_color() << std::endl;
        throw ArraySizeError(msg.str());
    }
}

// allocate GPU resources, this function will be executed by each fork
int
initCudaConstandTex(
        SHM_DATA     *datavols,        // host-side data arrays
        MONO_KERNELS *mono,            // spectrum data
        CONSTANTS    *constants,       // calculation paramaters/information
        int          nrays,            // total number of rays
        int          deviceid,         // non-CUDA device index
        int          gpuid,            // CUDA GPU ID assigned to this thread
        bool         verbose=false,    // output switch
        bool         timing=false,     // output switch
        bool         debug_write=false // write debug data to file
) {
    CudaTimer timer;
    if (timing) { timer.start(); }
    if (verbose) {
        cudaMemInfo meminfo = query_device_memory();
        meminfo.gpuid = gpuid;
            meminfo.print_available();
    }

    ///////////////// CONVOLUTION KERNEL / BEAM SPECTRUM / RADIAL BOUNDARY ////////////////////////////
    cudaChannelFormatDesc floatChannelDesc = cudaCreateChannelDesc<float>();

    // KERNEL WEIGHTS: copy to device memory cudaArray and bind texture reference
    checkCudaErrors( cudaMallocArray( &d_kern_array[deviceid], &floatChannelDesc, constants->nradii, constants->ntheta ) );
    checkCudaErrors( cudaMemcpyToArray(d_kern_array[deviceid], 0, 0, datavols->kernel_array, constants->nradii * constants->ntheta * sizeof(float), cudaMemcpyHostToDevice));
    makeTexObject<cudaArray>(&texKern[deviceid], d_kern_array[deviceid], 2, cudaAddressModeClamp, cudaFilterModeLinear);

    // spectrum data needed for TERMA calculation & normalization
    float *spectrum;
    spectrum = new float[4*mono->nkernels];
    float kermac0 = 0.0f, terma0 = 0.0f;
    for (int e=0; e<mono->nkernels; e++)
    {
        // calculate T and Kc at zero depth for use in beam hardening correction for terma
        // see notes in Siddon implementation
        // zero-depth terma/kerma used in correction of terma for polyenergetic kernel at depth
        kermac0 += mono->fluence[e]*mono->energy[e]*mono->mu_en[e];
        terma0 += mono->fluence[e]*mono->energy[e]*mono->mu[e];

        spectrum[e                   ] = mono->fluence[e];
        spectrum[e +   mono->nkernels] = mono->energy[e];
        spectrum[e + 2*mono->nkernels] = mono->mu_en[e];
        spectrum[e + 3*mono->nkernels] = mono->mu[e];
    }
    /* constants->beamhard_correct = terma0/kermac0; */
    constants->beamhard_correct = 1.f;

    // bind spectrum data to texture memory (with nearest-style fetching - no interpolation)
    // allocated as (nkernels x 4) dim matrix where:
    //   -- (:, 1): fluence
    //   -- (:, 2): energy
    //   -- (:, 3): mu_en
    //   -- (:, 4): mu
    // used in cudaSiddon for terma volume calculation
    checkCudaErrors( cudaMallocArray( &d_spectrum_array[deviceid], &floatChannelDesc, mono->nkernels, 4 ) );
    checkCudaErrors( cudaMemcpyToArray(d_spectrum_array[deviceid], 0, 0, spectrum, 4 * mono->nkernels * sizeof(float), cudaMemcpyHostToDevice));
    makeTexObject<cudaArray>(&texSpectrum[deviceid], d_spectrum_array[deviceid], 2, cudaAddressModeClamp, cudaFilterModePoint);

    // copy mean radii of kernel data to constant memory (first convert radial bounds to mean radii)
    /* float* mean_radii = new float[constants->nradii]; */
    /* memset(mean_radii, 0, constants->nradii*sizeof(float)); */
    /* mean_radii[0] = 0.5*datavols->radial_boundary[0]; */
    /* for (int rr=1; rr<constants->nradii; rr++) { */
    /*     mean_radii[rr] = 0.5*(datavols->radial_boundary[rr]+datavols->radial_boundary[rr-1]); */
    /* } */
    checkCudaErrors( cudaMemcpyToSymbol( KERN_RADII, datavols->radial_boundary, constants->nradii*sizeof(float) ) );
    /* delete[] mean_radii; */
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////// Bind Density to 3D Texture Array //////////////////////////////////////////////
    cudaExtent volExtent = make_cudaExtent(
            constants->size.x,
            constants->size.y,
            constants->size.z
            );
    cudaMalloc3DArray(&d_dens_array[deviceid], &floatChannelDesc, volExtent);
    cudaMemcpy3DParms CopyParams = {0};
    CopyParams.srcPtr   = make_cudaPitchedPtr(datavols->density, volExtent.width*sizeof(float), volExtent.width, volExtent.height);
    CopyParams.dstArray = d_dens_array[deviceid];
    CopyParams.extent   = volExtent;
    CopyParams.kind     = cudaMemcpyHostToDevice;
    cudaMemcpy3DAsync(&CopyParams);

    makeTexObject<cudaArray>(&texDens[deviceid], d_dens_array[deviceid], 3, cudaAddressModeBorder, cudaFilterModeLinear);
    /////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////// Find extent of BEV data volumes //////////////////////////////////////////////
    // Currently this is set to an isotropic box with side length == maximum diagonal length of calc_bbox
    // This will depend on the rotation angles for beam and convolution rays
    {
        // TODO: Do better job of estimating this
        float min_vox_size = min(constants->rev_latspacing, constants->rev_longspacing);
        unsigned int max_diag_len = ceilf(length(make_float3(constants->calc_bbox_size))*constants->voxel.x/min_vox_size);
        unsigned int buffer = 100; // sometimes rev->size exceeds estimate here by ~20
        constants->max_rev_size = make_uint3(
                max_diag_len + 2*(REV_PAD/min_vox_size + constants->penumbra/min_vox_size) + buffer,
                max_diag_len + 2*(REV_PAD/min_vox_size + constants->penumbra/min_vox_size) + buffer,
                max_diag_len + 2*(REV_PAD/min_vox_size + constants->penumbra/min_vox_size) + buffer
                );
    }
    printf("  Full BEV Volume Dimensions: %d x %d x %d\n",constants->max_rev_size.x,constants->max_rev_size.y,constants->max_rev_size.z);

    // Following memory operations are allocated on their assigned GPU

    // allocate intermediate data volume according to BEV extents // memsets occur later, just before use of each in kernel looping
    int revSize = constants->max_rev_size.x*constants->max_rev_size.y*constants->max_rev_size.z*sizeof(float);
    checkCudaErrors( cudaMalloc( (void**)&device_data[deviceid].bevDens,  revSize ) );
    checkCudaErrors( cudaMalloc( (void**)&device_data[deviceid].bevTerma, revSize ) );


    // Cuda Array for terma texture object fetching
    cudaExtent calc_bbox_extent = make_cudaExtent(constants->calc_bbox_size.x, constants->calc_bbox_size.y, constants->calc_bbox_size.z);
    cudaMalloc3DArray(&device_data[deviceid].term_Array, &floatChannelDesc, calc_bbox_extent);

    // allocate 3D CUDA Array to find with surface and texture objects
    cudaExtent doseExtent = make_cudaExtent( constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
    cudaMalloc3DArray(&(device_data[deviceid].dose_Array), &floatChannelDesc, doseExtent, cudaArraySurfaceLoadStore);

    // generate texture/surface object for reading/writing terma and dose data in kernels
    makeTexObject<cudaArray>(&(device_data[deviceid].texTerma), device_data[deviceid].term_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);
    makeTexObject<cudaArray>(&(device_data[deviceid].texDose),  device_data[deviceid].dose_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);
    makeSurfObject(&(device_data[deviceid].surfDose), device_data[deviceid].dose_Array);
    //////////////////////////////////////////////////////////////////////////////////////////////

    if (verbose) {
        cudaMemInfo meminfo = query_device_memory();
        meminfo.gpuid = gpuid;
            meminfo.print_available();
    }
    if (timing) {
        timer.stop_print_time_elapsed("GPU mem allocation");
    }

    return 1;
}

// performs the convolution for each beam
int
radconvolveTexture(
        MONO_KERNELS     *mono,
        CONSTANTS        *constants,
        BEAM             *device_beam_arr,
        int              nbeams,
        int              nrays,
        int              deviceid,      // non-CUDA device index (always start from 0, increment by 1)
        int              gpuid,         // CUDA device index (could skip ints for incompatible devices)
        bool             verbose=false,
        bool             timing=false,
        bool             debugwrite=false
) {
    cudaDeviceProp devProp;
    checkCudaErrors(cudaGetDeviceProperties(&devProp, gpuid));

    float3 rev_voxelsize = {constants->rev_latspacing, constants->rev_longspacing, constants->rev_latspacing};

    // set up arrays for dynamic GPU resource allocation
    // dependent on dimensions of BEV data per convolution ray
    dim3 *rayGrid = new dim3[nrays];
    dim3 *conBlock = new dim3[nrays];
    dim3 *conGrid = new dim3[nrays];
    unsigned int *memsize = new unsigned int[nrays];

    // Calculate cuda execution block/grid sizes
    unsigned int dataSize = constants->nvoxels();
    unsigned int calcDataSize = constants->bbox_nvoxels();
    unsigned int revSize = constants->max_rev_size.x * constants->max_rev_size.y * constants->max_rev_size.z;

    dim3 siddonBlock;
    dim3 siddonGrid;
    if (constants->ss_factor>1) {
        uint nsubvoxels = powf(constants->ss_factor, 3);
        uint siddonThreadLimit = 256;
        uint tilesize = floorf(sqrtf(siddonThreadLimit/nsubvoxels));
        siddonBlock = dim3{nsubvoxels, tilesize, tilesize};
    } else {
        siddonBlock = dim3{1, 32, 3};
    }
    siddonGrid = dim3(
            constants->calc_bbox_size.z,
            static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.x)/siddonBlock.y)),
            static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.y)/siddonBlock.z))
            );
    int siddonSharedMem = siddonBlock.x*siddonBlock.y*siddonBlock.z*sizeof(float);
    if (verbose && extra_verbose && deviceid == 0) {
        printf("Terma Supersampling: %s\n", (constants->ss_factor>1)?"enabled":"disabled");
        printf(" ### siddonGrid:  %3d x %3d x %3d  |  siddonBlock:  %3d x %3d x %3d | Shared: %d bytes\n",
                siddonGrid.x, siddonGrid.y, siddonGrid.z,
                siddonBlock.x, siddonBlock.y, siddonBlock.z,
                siddonSharedMem
              );
    }

    // create thread stream(s)
    cudaStream_t beamstream;
    cudaStreamCreate(&beamstream);

    CudaTimer timer(beamstream);
    if (timing) { timer.start(); }

    /////////////// Declare Reusable Memory ///////////////////
    // allocate terma and dose data volumes // memsets performed at start of each beam computation
    checkCudaErrors( cudaMalloc( (void**) &device_data[deviceid].terma, calcDataSize * sizeof(float) ) );
    checkCudaErrors( cudaMalloc( (void**) &device_data[deviceid].dose, dataSize * sizeof(float) ) );

    float *d_fluence_map; // device storage of beamlet fluence weights/intensities
    checkCudaErrors( cudaMalloc( (void**) &d_fluence_map, device_beam_arr[0].fmap_size.x*device_beam_arr[0].fmap_size.y*sizeof(float) ));

    float *h_terma;  // host (pinned) memory for finding cropped terma volume and termaExtent/lim_min/lim_max
    checkCudaErrors( cudaHostAlloc( (void **) &h_terma, calcDataSize * sizeof(float), 0) );

    float *h_beam_dose;  // host (pinned) memory for writing single beam dose volume to file
    checkCudaErrors( cudaHostAlloc( (void **) &h_beam_dose, dataSize * sizeof(float), 0) );
    ///////////////////////////////////////////////////////////

    if (timing) { timer.reset_print_time_elapsed("Convolution setup"); }

    // for each beam assigned to this GPU
    for (int dc=0; dc<nbeams; dc++) {
        BEAM* this_beam = &device_beam_arr[dc];
        if (verbose) {
            std::cout<<"| device#"<<gpuid<<" beam#"<<dc<<" deviceid#"<<deviceid<<" | beam="<<*this_beam<<std::endl;
        }
        if (timing) { timer.start(); }

        // collect per-beam memsets
        // TODO: some memsets may be unnecessary if memcpy fills array with new vals anyways
        checkCudaErrors( cudaMemsetAsync(device_data[deviceid].dose,        0, dataSize * sizeof(float), beamstream) );
        checkCudaErrors( cudaMemsetAsync(device_data[deviceid].terma,       0, calcDataSize * sizeof(float), beamstream) );
        checkCudaErrors( cudaMemsetAsync(h_terma,                            0, calcDataSize * sizeof(float), beamstream) );
        // copy binary fluence mask to GPU
        checkCudaErrors( cudaMemcpyAsync(d_fluence_map, this_beam->fluence_map.data(), this_beam->fmap_size.x*this_beam->fmap_size.y*sizeof(float), cudaMemcpyHostToDevice, beamstream) );

        // calculate terma
        cudaSiddon<<< siddonGrid, siddonBlock, siddonSharedMem, beamstream >>>(
                device_data[deviceid].terma,
                d_fluence_map,
                constants->start,
                make_float3(constants->size),
                constants->voxel,
                constants->calc_bbox_start,
                constants->calc_bbox_size,
                constants->beamhard_correct,
                this_beam->source,
                this_beam->direction,
                this_beam->isocenter,
                this_beam->sad,
                this_beam->beamlet_size,
                this_beam->fmap_size,
                this_beam->azimuth,
                this_beam->zenith,
                this_beam->coll,
                mono->nkernels,
                texDens[deviceid],
                texSpectrum[deviceid]
                );
        checkCudaErrors( cudaStreamSynchronize(beamstream) ); // wait for cudaSiddon to complete
        getLastCudaError("Kernel execution failed");

        // find terma extent
        checkCudaErrors( cudaMemcpyAsync( h_terma, device_data[deviceid].terma, calcDataSize * sizeof(float), cudaMemcpyDeviceToHost, beamstream ) );
        checkCudaErrors( cudaStreamSynchronize(beamstream) ); // wait for cudaSiddon to complete

        // write debug output (terma) to file
        if (debugwrite) {
            // save terma to file
            char checkchar[50];
            sprintf(checkchar,"terma-d%d-b%d", gpuid, dc);
            write_debug_data( h_terma, constants->calc_bbox_size, checkchar, verbose);
        }

        // find coords of terma bounding box start/end in RCS
        float3 lim_min, lim_max;
        uint3 bev_cropbox_size;
        { // local scope
            lim_min = lim_max = make_float3(0.f,0.f,0.f);
            // temporary storage of BEV indices for min and max voxel (rotated to BEV orientation rather than RCS)
            float3 b_min_indices, b_max_indices;
            b_max_indices = make_float3(std::numeric_limits<float>::min());
            b_min_indices = make_float3(std::numeric_limits<float>::max());

            // TODO: consider implementing GPU kernel (reduction) to avoid device-to-host memcpy latencies (~22ms on CPU; ~4% of beam exec time)
            // step through Global coord sys and if non-zero check bounds in BEV coord sys
            // report limits in RCS coord sys
            for (int z=0; z<constants->calc_bbox_size.z; z++) {
                for (int y=0; y<constants->calc_bbox_size.y; y++) {
                    for (int x=0; x<constants->calc_bbox_size.x; x++) {
                        int vid = x + constants->calc_bbox_size.x * (y + constants->calc_bbox_size.y * z);
                        if (h_terma[vid] != 0) {
                            // rotate into BEV coord sys
                            float3 b_indices = rotateBeamAtOriginRHS(make_float3(x,y,z), this_beam->azimuth, this_beam->zenith, this_beam->coll);
                            /* std::cout << FORMAT_3VEC((int3{x,y,z}))<<" --> " <<FORMAT_3VEC(b_indices) << std::endl; */
                            b_min_indices.x = fminf(b_min_indices.x, b_indices.x);
                            b_min_indices.y = fminf(b_min_indices.y, b_indices.y);
                            b_min_indices.z = fminf(b_min_indices.z, b_indices.z);
                            b_max_indices.x = fmaxf(b_max_indices.x, b_indices.x);
                            b_max_indices.y = fmaxf(b_max_indices.y, b_indices.y);
                            b_max_indices.z = fmaxf(b_max_indices.z, b_indices.z);
                        }
                    }
                }
            }
            // increase terma extent by penumbra size in BEV coord sys
            b_min_indices.x -= (constants->penumbra / constants->rev_latspacing);
            b_min_indices.y -= (constants->penumbra / constants->rev_longspacing);
            b_min_indices.z -= (constants->penumbra / constants->rev_latspacing);
            b_max_indices.x += (constants->penumbra / constants->rev_latspacing);
            b_max_indices.y += (constants->penumbra / constants->rev_longspacing);
            b_max_indices.z += (constants->penumbra / constants->rev_latspacing);

            // size of cropped bev volume aligned to BEV coordinate system (beam orientation in +Y)
            bev_cropbox_size = make_uint3(
                    ceilf(abs(b_max_indices.x - b_min_indices.x)),
                    ceilf(abs(b_max_indices.y - b_min_indices.y)),
                    ceilf(abs(b_max_indices.z - b_min_indices.z))
                    );

            // convert bev extents to RCS coords where bounding box is aligned with BEV coordinate system
            float3 temp_min = inverseRotateBeamAtOriginRHS(b_min_indices, this_beam->azimuth, this_beam->zenith, this_beam->coll);
            float3 temp_max = inverseRotateBeamAtOriginRHS(b_max_indices, this_beam->azimuth, this_beam->zenith, this_beam->coll);
            lim_min = make_float3(
                    constants->start.x + (temp_min.x+constants->calc_bbox_start.x)*constants->voxel.x,
                    constants->start.y + (temp_min.y+constants->calc_bbox_start.y)*constants->voxel.y,
                    constants->start.z + (temp_min.z+constants->calc_bbox_start.z)*constants->voxel.z
                    );
            lim_max = make_float3(
                    constants->start.x + (temp_max.x+constants->calc_bbox_start.x)*constants->voxel.x,
                    constants->start.y + (temp_max.y+constants->calc_bbox_start.y)*constants->voxel.y,
                    constants->start.z + (temp_max.z+constants->calc_bbox_start.z)*constants->voxel.z
                    );

            if (verbose) {
                printf("Cropbox Calculation limits (scanner coords):  X:%0.2f-->%0.2f, Y:%0.2f-->%0.2f, Z:%0.2f-->%0.2f\n", lim_min.x, lim_max.x, lim_min.y, lim_max.y, lim_min.z, lim_max.z);
                printf("Cropbox Dimensions (BEV-aligned): %d x %d x %d\n", bev_cropbox_size.x, bev_cropbox_size.y, bev_cropbox_size.z );
                /* std::cout << "b_min_indices: "<<FORMAT_3VEC(b_min_indices)<<" b_max_indices: "<<FORMAT_3VEC(b_max_indices) << std::endl; */
            }
            if (bev_cropbox_size.x <= 0 || bev_cropbox_size.x >= 1000.f ||
                bev_cropbox_size.y <= 0 || bev_cropbox_size.y >= 1000.f ||
                bev_cropbox_size.z <= 0 || bev_cropbox_size.z >= 1000.f) {
                throw std::runtime_error("bev_cropbox_size error\n");
            }
        }

        // copy terma volume to 3D cuda array, and bind to texture with local scope
        // // term_array[deviceid] allocated outside of beam loop
        { // local scope
            cudaExtent calc_bbox_extent = make_cudaExtent(constants->calc_bbox_size.x, constants->calc_bbox_size.y, constants->calc_bbox_size.z);
            cudaMemcpy3DParms CopyParams = {0};
            CopyParams.srcPtr   = make_cudaPitchedPtr((void*)device_data[deviceid].terma, calc_bbox_extent.width*sizeof(float), calc_bbox_extent.width, calc_bbox_extent.height);
            CopyParams.dstArray = device_data[deviceid].term_Array;
            CopyParams.extent   = calc_bbox_extent;
            CopyParams.kind     = cudaMemcpyDeviceToDevice;
            checkCudaErrors( cudaMemcpy3DAsync(&CopyParams, beamstream) );
        }


        // Prepare kernel launch parameters
        dim3 rayBlock(32, 3, 1);

        dim3 xyzBlock = dim3(TILE_DIM_X, TILE_DIM_Y, 1);
        dim3 xyzGrid = dim3(
                static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.x)/xyzBlock.x)),
                static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.y)/xyzBlock.y)),
                static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.z)/xyzBlock.z))
                );
        if (verbose && extra_verbose) {
            printf(" ### xyzGrid:   %3d x %3d x %3d  |  xyzBlock:   %3d x %3d x %3d\n",
                    xyzGrid.x, xyzGrid.y, xyzGrid.z,
                    xyzBlock.x, xyzBlock.y, xyzBlock.z
                  );
        }

        // holds dimensions and data for the active volume when viewed at each convolution angle
        // used for dynamic GPU memory and resource allocation
        REV_DATA* rev = new REV_DATA[nrays];

        // find the dimensions of the transformed data in BEV/NVB coordinate system for all rays
        // rev->min_coords and bev->max_coords will be an expansion of lim_min and lim_max based on REV rotation
        for (int i=0; i<nrays; i++) {
            // compute angles from idx
            if (verbose && extra_verbose) { printf(" Ray #%d\n", i+1); }
            try {
                findREV(&rev[i],
                        constants,
                        this_beam,
                        make_float3(0.f),
                        make_float3(bev_cropbox_size)*rev_voxelsize,
                        constants->get_theta_from_index(i),
                        constants->get_phi_from_index(i),
                        (verbose && extra_verbose)
                        );
            } catch (const ArraySizeError& e) {
                std::cerr << e.what();
                if (!debugwrite) { throw e; }
            }

            // PAY ATTENTION: To keep memory accesses coalesced where possible coordinates and data layout axes are not the same
            // as such, kernel launch param axes also do not match. coordinate axes XYZ are mapped to device memory layout axes ZXY
            // determine block and grid sizes for all convolution functions for all convolution rays
            rayGrid[i] = make_uint3(
                    static_cast<unsigned int>(ceilf(static_cast<float>(rev[i].size.x)/rayBlock.x)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(rev[i].size.y)/rayBlock.y)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(rev[i].size.z)/rayBlock.z))
                    );

            // x    # samples along conv. ray direction
            // y|z: # parallel rays in REV
            // TODO: be careful when rev[rr].size.x approaches #threads-per-block limit (1024 for CC3.0+)
            //       In most cases this should be below 300-400 though
            if (rev[i].size.x > 1024) { throw std::runtime_error("size of rev (x-dim) exceeds maximum #threads allowed. Contact dev to resolve."); }
            int overwrite = 2; // buffer beyond used rev block, where convolve kernel will also write zeros before proceeding
            conBlock[i].x = rev[i].size.x + overwrite;
            conGrid[i].y  = rev[i].size.y + overwrite;
            conGrid[i].z  = rev[i].size.z + overwrite;
            conBlock[i].y = conBlock[i].z = conGrid[i].x = 1;
            // dyn. alloc. shared mem size
            memsize[i] = 2 * (conBlock[i].x - overwrite) * sizeof(float);


            if (verbose && extra_verbose) {
                printf(" ### rayGrid:   %3d x %3d x %3d  |  rayBlock:   %3d x %3d x %3d\n",
                        rayGrid[i].x, rayGrid[i].y, rayGrid[i].z,
                        rayBlock.x, rayBlock.y, rayBlock.z
                      );
                printf(" ### conGrid:   %3d x %3d x %3d  |  conBlock:   %3d x %3d x %3d\n\n",
                        conGrid[i].x, conGrid[i].y, conGrid[i].z,
                        conBlock[i].x, conBlock[i].y, conBlock[i].z
                      );
            }
        }

        if (timing) { timer.restart_print_time_elapsed("Beam setup"); }

        int debug_ray = constants->ntheta/2+1;
        for (int ray_idx=0; ray_idx<nrays; ray_idx++) {
            // get kernel rotation angles once instead of for each kernel thread launch
            int kern_wt_idx =  constants->get_kernel_theta_index(ray_idx);
            float kern_theta = constants->get_theta_from_index(ray_idx);
            float kern_phi =   constants->get_phi_from_index(ray_idx);

            // kernels are launched with # of threads according to pre-calculated dimensions of BEV data
            // for this beam and ray direction
            // results are stored on this threads GPU global mem

            // initialize resampling volumes
            checkCudaErrors( cudaMemsetAsync( device_data[deviceid].bevDens,      0, revSize * sizeof(float), beamstream ) );
            checkCudaErrors( cudaMemsetAsync( device_data[deviceid].bevTerma,     0, revSize * sizeof(float), beamstream) );

            // ray trace along kernel CCCS ray (for each direction) and resample to REV volume
            cudaRaytrace <<< rayGrid[ray_idx], rayBlock, 0, beamstream >>> (
                    device_data[deviceid].bevDens,  // BEV output
                    device_data[deviceid].bevTerma, // BEV output
                    this_beam->azimuth,
                    this_beam->zenith,
                    this_beam->coll,
                    kern_theta, kern_phi,                // ray direction
                    rev[ray_idx].min_coords,           // rev box start indices
                    lim_min,
                    rev[ray_idx].size,                   // REV volume dims
                    bev_cropbox_size,
                    constants->max_rev_size,
                    constants->start,
                    constants->voxel,
                    rev_voxelsize,
                    make_float3(constants->calc_bbox_start),
                    make_float3(constants->calc_bbox_start + constants->calc_bbox_size),
                    texDens[deviceid],                   // XYZ input
                    device_data[deviceid].texTerma      // XYZ input
                    );

            if (debugwrite && extra_debug && ray_idx <= debug_ray) {
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
                float* out_array = new float[revSize * sizeof(float) ];
                memset((void*)out_array, 0, revSize * sizeof(float));
                checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].bevDens, revSize*sizeof(float), cudaMemcpyDeviceToHost));

                char checkchar[40];
                sprintf(checkchar,"debug_dens-d%d-b%d-r%d", gpuid, dc, ray_idx);
                write_debug_data(out_array, constants->max_rev_size, checkchar, verbose);
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
                /////////////////////////////////////////////////////
                memset((void*)out_array, 0, revSize * sizeof(float));
                checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].bevTerma, revSize*sizeof(float), cudaMemcpyDeviceToHost));

                sprintf(checkchar,"debug_terma-d%d-b%d-r%d", gpuid, dc, ray_idx);
                write_debug_data(out_array, constants->max_rev_size, checkchar, verbose);
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
            }

            // perform dose calcluation (CCCS) w/ heterogeneity correction in BEV volume
            RowConvolve <<< conGrid[ray_idx], conBlock[ray_idx], memsize[ray_idx], beamstream >>> (
                    device_data[deviceid].bevDens,  // BEV-trans input
                    device_data[deviceid].bevTerma, // BEV-trans input
                    device_data[deviceid].surfDose, // BEV dose output
                    (float)kern_wt_idx,              // kernel theta index (for sampling texKern)
                    rev[ray_idx].size,               // BEV volume dims
                    constants->max_rev_size,
                    constants->rev_longspacing,
                    constants->nradii,
                    constants->ntheta,
                    constants->nphi,
                    constants->penumbra,
                    texKern[deviceid]                // kernel weights
                    );

            if (debugwrite && extra_debug && ray_idx <= debug_ray) {
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
                float* out_array = new float[revSize * sizeof(float) ];
                memset((void*)out_array, 0, revSize * sizeof(float));
                cudaExtent doseExtent = make_cudaExtent( constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
                cudaMemcpy3DParms CopyParams = {0};
                CopyParams.srcArray = device_data[deviceid].dose_Array;
                CopyParams.dstPtr   = make_cudaPitchedPtr(out_array, sizeof(float)*doseExtent.width, doseExtent.width, doseExtent.height);
                CopyParams.extent   = doseExtent;
                CopyParams.kind     = cudaMemcpyDeviceToHost;
                checkCudaErrors(cudaMemcpy3D(&CopyParams));

                char checkchar[40];
                sprintf(checkchar,"debug3_dose-d%d-b%d-r%d", gpuid, dc, ray_idx);
                write_debug_data(out_array, constants->max_rev_size, checkchar, verbose);
                checkCudaErrors(cudaStreamSynchronize(beamstream));
            }

            // resample REV dose back to dicom coordinate system (XYZ)
            // dose for this ray (direction) is stored in device_data[deviceid].dose
            // all dose volumes calculated for this device are automatically summed together here
            float3 revStart_inRCS = lim_min + inverseRotateBeamAtOriginRHS(rev[ray_idx].min_coords, this_beam->azimuth, this_beam->zenith, this_beam->coll);
            REV2XYZdose <<< xyzGrid, xyzBlock, 0, beamstream >>> (
                    device_data[deviceid].dose,    // XYZ dose output
                    device_data[deviceid].texDose, // BEV dose input
                    this_beam->azimuth,
                    this_beam->zenith,
                    this_beam->coll,
                    kern_theta, kern_phi,           // ray direction
                    constants->start,
                    constants->voxel,
                    rev_voxelsize,
                    constants->size,
                    constants->calc_bbox_start,
                    constants->calc_bbox_start+constants->calc_bbox_size,
                    rev[ray_idx].min_coords,
                    rev[ray_idx].size,
                    lim_min
                    );

        }
        getLastCudaError("Kernel execution failed");
        checkCudaErrors( cudaStreamSynchronize(beamstream) );  // wait for all rays to be completed

        if (timing) {
            timer.stop();
            float cudatime = timer.time_elapsed();
            printf("-- %d.%d - Kernel Convolve Time: %4.3f msec (avg time: %4.3f msec per ray (of %d rays) --\n", gpuid, dc, cudatime, cudatime/nrays, nrays);
            printf("\n"); // blankline between beams of same device
        }

        // write dose to file
        if (debugwrite) {
            checkCudaErrors( cudaMemcpyAsync(h_beam_dose, device_data[deviceid].dose, dataSize * sizeof(float), cudaMemcpyDeviceToHost, beamstream));
            checkCudaErrors( cudaStreamSynchronize(beamstream) );
            std::ostringstream temp_stream;
            temp_stream<<"dose-d"<<gpuid<<"-b"<<dc;
            std::string dose_outfile = temp_stream.str();
            write_debug_data(h_beam_dose, constants->size, dose_outfile.c_str(), verbose);
        }

        delete [] rev;
    }
    /* checkCudaErrors( cudaStreamSynchronize(beamstream) );  // wait for all beams to be completed */
    if (timing) { timer.start(); }

    // cleanup thread memory
    delete [] rayGrid;
    delete [] conBlock;
    delete [] conGrid;
    delete [] memsize;
    checkCudaErrors( cudaFree(d_fluence_map) );
    checkCudaErrors( cudaFreeHost(h_terma) );
    checkCudaErrors( cudaFreeHost(h_beam_dose) );
    checkCudaErrors( cudaFreeArray(device_data[deviceid].term_Array) );

    if (timing) { timer.stop_print_time_elapsed("Convolution cleanup"); }

    return 1;
}

// free GPU resources for each fork
int
freeCudaTexture(int ndevices, int *gpuid_arr)
{
    for (int deviceid=0; deviceid<ndevices; deviceid++) {
        int gpuid = gpuid_arr[deviceid];
        checkCudaErrors( cudaSetDevice(gpuid) );
        checkCudaErrors(cudaFree(device_data[deviceid].bevDens));
        checkCudaErrors(cudaFree(device_data[deviceid].bevTerma));
        checkCudaErrors(cudaDestroySurfaceObject(device_data[deviceid].surfDose));
        checkCudaErrors(cudaDestroyTextureObject(device_data[deviceid].texDose));
        checkCudaErrors(cudaDestroyTextureObject(device_data[deviceid].texTerma));
        checkCudaErrors(cudaFreeArray(device_data[deviceid].dose_Array));
        checkCudaErrors(cudaFree(device_data[deviceid].dose));
        checkCudaErrors(cudaFree(device_data[deviceid].terma));

        checkCudaErrors(cudaDestroyTextureObject(texDens[deviceid]));
        checkCudaErrors(cudaFreeArray(d_dens_array[deviceid]));

        checkCudaErrors(cudaDestroyTextureObject(texSpectrum[deviceid]));
        checkCudaErrors(cudaFreeArray(d_spectrum_array[deviceid]));

        checkCudaErrors(cudaDestroyTextureObject(texKern[deviceid]));
        checkCudaErrors(cudaFreeArray(d_kern_array[deviceid]));
    }

    return 1;
}

// sum dose results from each GPU, only run by one thread
int
cuda_sum_device_float_dose(
        float *dose,
        int ndevices,
        CONSTANTS *constants,
        bool timing,
        bool verbose
) {
    CudaTimer timer_task, timer_total;
    if (timing) {
        timer_task.start();
        timer_total.start();
    }

    dim3 sumBlock(_MAX_THREADS_PER_BLOCK_/2);
    //unsigned int convo_size = constants->termaSize.x * constants->termaSize.y * constants->termaSize.z;
    unsigned int convo_size = constants->size.x * constants->size.y * constants->size.z;

    unsigned int xdim = convo_size / sumBlock.x;
    if (xdim < 1) { xdim = 1; }
    unsigned int ydim = 1;
    if (convo_size % sumBlock.x > 0) { xdim++; }
    if (xdim >= _MAX_GRID_SIZE_) {
        ydim = xdim / _MAX_GRID_SIZE_;
        if (xdim % _MAX_GRID_SIZE_ > 0) { ydim++; }
        xdim = _MAX_GRID_SIZE_;
    }
    dim3 sumGrid(xdim,ydim);

    unsigned int data_size = constants->size.x * constants->size.y * constants->size.z;
    float *d_dose;
    checkCudaErrors( cudaMalloc( (void**) &d_dose, data_size * sizeof(float) ) );
    checkCudaErrors( cudaMemset( d_dose, 0, data_size * sizeof(float) ) );

    if (timing) { timer_task.restart_print_time_elapsed("Dose summation setup"); }

    // for each additional GPU used by the server, the dose is added to the parent's GPU-0 results
    for (int deviceid=0; deviceid<ndevices; deviceid++) {
        cudaSum <<< sumGrid, sumBlock >>> ( d_dose, device_data[deviceid].dose, data_size );
    }

    /* // use thrust library to find max */
    /* // map GPU dose volume to thrust pointer */
    /* thrust::device_ptr<float> dvec( d_dose ); */

    /* // find maximum value of dose */
    /* // would normally be used to normalize dose results */
    /* // but for mgcs framework, normalization is performed once results have been accumulated from all machines */
    /* float maxval = thrust::reduce( dvec, dvec + data_size, (float)0, thrust::maximum<float>() ); */
    /* if (verbose) */
    /*     printf(" Maximum Summed Dose: %f\n", maxval ); */

    // copy summed result back to host
    // dose to be normalized by host
    checkCudaErrors( cudaMemcpy( dose, d_dose, data_size * sizeof(float), cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaFree ( d_dose ) );

    if (timing) {
        timer_task.stop_print_time_elapsed("Dose sum execution");
        timer_total.stop_print_time_elapsed("Total dose summation");
    }

    return 1;
}
