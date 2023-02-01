#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring> // memset
#include <memory>
#include <limits>

#include "profile.h"
#include "DoseCalcAlgo/cudaSiddon.cuh"
#include "nvbbRayConvolve_device.cu" // device memory pointers & kernel functions
#include "dosecalc_defs.h"  // MAXIMUM_DEVICE_COUNT
#include "server/brain_defs.h" // CONSTANTS, SHM_DATA ...

#include "Utilities/multithreading.h"
#include "Utilities/timing.h"
#include "CudaUtilities/memory_manager.cuh"
#include "CudaUtilities/cuda_timing.cuh"
#include "CudaUtilities/manage_gpu.cuh"
#include "CudaUtilities/make_tex_surf.cuh"
#include "CudaUtilities/geometry.cuh" // coordinate rotations
#include "CudaUtilities/array_sum.cuh"
#include "DoseCalcIO/dosecalcio.h"
#include "sparsify_manager.h"

// Sparsify Worker Queue from main.cu shared across threads
extern SparsifyManager   sparsifymanager;
extern ROIMaskList       roi_list;
extern CONSTANTS*        constants;

extern bool extra_verbose; // print findREV verbose output
extern bool extra_debug;    // print between-kernel-execution data structures (dose, terma, density)
extern bool TERMA_ONLY;

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
    /* float *terma; */
    float *dose;                    // storage of final dose for each beam
    float *device_dose;             // cleared and set for each beam calculation on each device
    float *bevDens;
    float *bevTerma;
    cudaArray *term_Array;
    cudaArray *dose_Array;
    cudaTextureObject_t texTerma;   // texture object created/destroyed for each beam - array mem is once initialized for all beams
    cudaTextureObject_t texDose;
    cudaSurfaceObject_t surfDose;
};


// Unique for each device
DEVICE_CONV_DATA  device_data        [ MAXIMUM_DEVICE_COUNT ] = {};       // holds cuda execution device specific data
float*            device_dose_accum  [ MAXIMUM_DEVICE_COUNT ] = {};       // holds accumulated dose for each device (summed in main.cu later)
MemoryManager     memory_managers    [ MAXIMUM_DEVICE_COUNT ] = {};

static std::vector<int> allocate_batch_sizes(int numbeamlets, int numbatches) {
    int min_batch_size = floor(numbeamlets/numbatches); // number of batches with size1
    // evenly allocate beamlets to batches
    std::vector<int> batch_sizes(numbatches, min_batch_size);
    int remaining = numbeamlets - min_batch_size*numbatches;
    int _mark = 0;
    while (remaining>0) {
        batch_sizes[(_mark++)%numbatches] += 1;
        remaining--;
    }
    return batch_sizes;
}

/* Get RCS coordinates of start/end voxels in REV coordinate system rotated to REV orientation
 * This is critical for converting between REV indices to RCS coords then RCS indices to sample terma/density in RCS
 */
static void update_extents(float3& currMin, float3& currMax, const float3& thisPt) {
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
        float theta;
        float phi;
        uint3 oldsize;
        REV_DATA* rev;
        ArraySizeError(float theta, float phi, uint3 oldsize, REV_DATA* rev) :
            std::runtime_error{""}, theta{theta}, phi{phi}, oldsize{oldsize}, rev{rev} {}
};
// calculate the extent of packedBEV in REV
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
    float3 b_diff = lim_max-lim_min;
    float3 pbev_min = make_float3(std::numeric_limits<float>::max());
    float3 pbev_max = make_float3(std::numeric_limits<float>::min());

    // for 8 corners of the data volume defined by RCS coordinates of BEV box, find position in REV coord sys
    for (int xoff=0; xoff<2; xoff++) {
        for (int yoff=0; yoff<2; yoff++) {
            for (int zoff=0; zoff<2; zoff++) {
                float3 input = make_float3(b_diff.x*xoff, b_diff.y*yoff, b_diff.z*zoff);
                // apply rotation for the convolution ray and evaluate extents in REV space
                float3 output = rotateKernelAtOriginRHS(input, theta, phi);
                // set beam's eye view extents
                update_extents(pbev_min, pbev_max, output);
    } } }

    // rotate limit coords into RCS space and shift relative to RCS origin
    // remember: Pillar_grid to REV orientation is XYZ -> ZXY and rev size reflects REV orientation
    rev->size = make_uint3(
            static_cast<unsigned int>( ceil(fabsf(pbev_max.y-pbev_min.y) / constants->rev_longspacing) ),
            static_cast<unsigned int>( ceil(fabsf(pbev_max.z-pbev_min.z) / constants->rev_latspacing) ),
            static_cast<unsigned int>( ceil(fabsf(pbev_max.x-pbev_min.x) / constants->rev_latspacing) ) );

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
    if (was_shrunk) { pbev_max -= make_float3(adjust.z, adjust.x, adjust.y); }

    // store limits in REV coordinates relative to rotated lim_min, lim_max
    rev->min_coords = pbev_min;
    rev->max_coords = pbev_max;

    if (verbose) {
        printf(" Theta               :  %5.1f deg\n"            , theta*180/PI );
        printf(" Phi                 :  %5.1f deg\n"            , phi*180/PI );
        printf(" REV->size (PG:YZX)  :    %5d x   %5d x   %5d\n", rev->size.x, rev->size.y, rev->size.z);
        printf(" pbev_min            :  %7.2f x %7.2f x %7.2f\n", rev->min_coords.x, rev->min_coords.y, rev->min_coords.z);
        printf(" pbev_max            :  %7.2f x %7.2f x %7.2f\n", rev->max_coords.x, rev->max_coords.y, rev->max_coords.z);
    }
    if (was_shrunk) {
        throw ArraySizeError(theta, phi, oldsize, rev);
    }
}

/* each beamlet is defined by the anchor coords corresponding to the points along the beamlet
 * central axis that are nearest and furthest from the beam source, and within the calc_bbox
 * This function returns the start/end anchor GCS coords given the beamlet index
 */
#define FCOMP(x, ii) *(((float*)&x)+ii) // return component of float vector by index
static void calcBeamletAnchors(
        float3&      start,            // out: start_anchor_coords
        float3&      end,              // out: end_anchor_coords
        float2&      beamletAngles,    // out: beam angles + beamlet divergence angles (.x: azi, .y: .zen)
        float3&      beamletIso,       // out: beamlet isocenter coords
        float3       src,              // beam src coords
        float3       iso,              // beam iso coords
        unsigned int beamlet_idx,      // fluence map linearized index
        float2       beamlet_size,     // [unit: cm]
        uint2        fmap_dims,        // number of beamlets along each axis
        float3       voxelsize,        // [unit: cm]
        float3       density_start,    // coords of start of density volume
        uint3        calc_bbox_start,  // nvoxel offset from density start
        uint3        calc_bbox_size,   // nvoxel size from calc_bbox_start
        float        gantry_angle_rad, // [unit: rad]
        float        couch_angle_rad,  // [unit: rad]
        float        coll_angle_rad,   // [unit: rad]
        int          verbose=false
) {
    // get beamlet central axis coords on fmap (beamlet-isocenter)
    unsigned int bx = beamlet_idx % fmap_dims.x;
    unsigned int by = beamlet_idx / fmap_dims.x;
    float3 bdiff = make_float3(
            beamlet_size.x * (-0.5f*fmap_dims.x + bx + 0.5f),
            0,
            beamlet_size.y * (-0.5f*fmap_dims.y + by + 0.5f) );
    beamletIso = iso + inverseRotateBeamAtOriginRHS(bdiff, gantry_angle_rad, couch_angle_rad, coll_angle_rad);
    // TODO: CHECK THIS (beamletangles)
    beamletAngles = make_float2(
            acosf(length(iso-src)/length(beamletIso-src)),
            atan2f(bdiff.z, bdiff.x) );

    // get coords of bbox limits
    float3 bbox_begin = density_start + voxelsize * make_float3(calc_bbox_start);
    float3 bbox_end   = bbox_begin + voxelsize * make_float3(calc_bbox_size);

    double distance_buffer = 0.05;

    // evaluate intersection with each of 6 calc_bbox faces
    float3 intersections[6] = {};
    start = make_float3(std::numeric_limits<float>::max());
    end = make_float3(std::numeric_limits<float>::min());
    double min_dist = std::numeric_limits<double>::max();
    double max_dist = std::numeric_limits<double>::min();
    for (int ii=0; ii<3; ii++) { // for x, y, z
        float3 diff_src2b = beamletIso - src;
        double denom = FCOMP(diff_src2b, ii);
        double alpha1 = ( FCOMP(bbox_begin,ii) - FCOMP(src,ii) ) / (denom+1e-12);
        double alpha2 = ( FCOMP(bbox_end,ii)   - FCOMP(src,ii) ) / (denom+1e-12);
        intersections[2*ii  ] = src + alpha1 * diff_src2b;
        intersections[2*ii+1] = src + alpha2 * diff_src2b;
        if (verbose>2) {
            for (int _jj=0; _jj<2; _jj++) {
                float3& _tmp = intersections[2*ii+_jj];
                std::cout << "bbox intersection #"<<(ii*2+_jj) <<": "<<_tmp.x << ", "<<_tmp.y<<", "<<_tmp.z<<")"<< std::endl;
            }
        }
    }

    // check for valid intersection with calc_bbox faces
    for (int ii=0; ii<3; ii++) { // for x, y, x
        for (int jj=0; jj<2; jj++) { // for each intersection (of 2 per axis)
            int idx = 2*ii+jj;
            // given intersection with one dimension, do other two dims occur within bounds of box?
            if (FCOMP(intersections[idx],(ii+1)%3)+distance_buffer >= FCOMP(bbox_begin,(ii+1)%3) && FCOMP(intersections[idx],(ii+1)%3)-distance_buffer <= FCOMP(bbox_end,(ii+1)%3) &&
                FCOMP(intersections[idx],(ii+2)%3)+distance_buffer >= FCOMP(bbox_begin,(ii+2)%3) && FCOMP(intersections[idx],(ii+2)%3)-distance_buffer <= FCOMP(bbox_end,(ii+2)%3) )
            {
                if (verbose>2) {
                    float3& _tmp = intersections[idx];
                    std::cout << "valid intersection found (#"<<idx<<"): "<<_tmp.x << ", "<<_tmp.y<<", "<<_tmp.z<<")" << std::endl;
                }
                float3 _inter = intersections[idx];
                if (length(_inter - src) < min_dist) {
                    start = _inter;
                    min_dist = length(start - src);
                }
                if (length(_inter - src) > max_dist) {
                    end = _inter;
                    max_dist = length(end - src);
                }
            }
        }
    }

    if (!(min_dist <std::numeric_limits<float>::max() && max_dist > std::numeric_limits<float>::min())) {
        std::ostringstream msg;
        msg << "Failed to determine beamlet anchor coordinates for beamlet #" + std::to_string(beamlet_idx);
        if (verbose > 1) {
        msg << std::endl <<
            "DEBUG INFO:" << std::endl <<
            "beamlet angles:  az:"<<beamletAngles.x*180.f/PI <<"; zen:"<<beamletAngles.y*180.f/PI<<" deg)" << std::endl <<
            "bdiff:      ("<<bdiff.x<<","<<bdiff.y<<","<<bdiff.z<<")" << std::endl <<
            "biso:       ("<<beamletIso.x<<","<<beamletIso.y<<","<<beamletIso.z<<")" << std::endl <<
            "start:      ("<<start.x<<","<<start.y<<","<<start.z<<")" << std::endl <<
            "end:        ("<<end.x<<","<<end.y<<","<<end.z<<")" << std::endl <<
            "bbox_begin: ("<<bbox_begin.x<<","<<bbox_begin.y<<","<<bbox_begin.z<<")" << std::endl <<
            "bbox_end:   ("<<bbox_end.x<<","<<bbox_end.y<<","<<bbox_end.z<<")" << std::endl;
        }
        throw std::runtime_error(msg.str());
    }
}

// copy constant and tex-ref data to device (executed by each host thread)
int initCudaConstandTex(
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
    MemoryManager& memmgr = memory_managers[deviceid];
    memmgr.setMemoryLimit(query_device_memory(gpuid).free*GB); // TODO: move this to main.cu later
    std::cout << "Memory Management (device: "<<gpuid<<") - Available device memory: "<<(memmgr.memoryRemaining()/GB)<<" GB"<<std::endl;

    // enlarge L1 cache size (shrinking shared mem for these kernel which is unused anyways)
    cudaFuncSetCacheConfig(cudaBeamletRaytrace, cudaFuncCache::cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(PackedREVtoBEVdose, cudaFuncCache::cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(UnpackBEVDosePillar, cudaFuncCache::cudaFuncCachePreferL1);

    ///////////////// CONVOLUTION KERNEL / BEAM SPECTRUM / RADIAL BOUNDARY ////////////////////////////
    cudaChannelFormatDesc floatChannelDesc = cudaCreateChannelDesc<float>();

    // KERNEL WEIGHTS: copy to device memory cudaArray and bind texture reference
    memmgr.cudaMallocArray( &d_kern_array[deviceid], &floatChannelDesc, constants->nradii, constants->ntheta );
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
    constants->beamhard_correct = terma0/kermac0;
    constants->beamhard_correct = 1.0f; // XXX

    // bind spectrum data to texture memory (with nearest-style fetching - no interpolation)
    // allocated as (nkernels x 4) dim matrix where:
    //   -- (:, 1): fluence
    //   -- (:, 2): energy
    //   -- (:, 3): mu_en
    //   -- (:, 4): mu
    // used in cudaSiddon for terma volume calculation
    memmgr.cudaMallocArray( &d_spectrum_array[deviceid], &floatChannelDesc, mono->nkernels, 4 );
    checkCudaErrors( cudaMemcpyToArray(d_spectrum_array[deviceid], 0, 0, spectrum, 4 * mono->nkernels * sizeof(float), cudaMemcpyHostToDevice));
    makeTexObject<cudaArray>(&texSpectrum[deviceid], d_spectrum_array[deviceid], 2, cudaAddressModeClamp, cudaFilterModePoint);
    delete [] spectrum;

    // copy mean radii of kernel data to constant memory (first convert radial bounds to mean radii)
    float* mean_radii = new float[constants->nradii];
    memset(mean_radii, 0, constants->nradii*sizeof(float));
    mean_radii[0] = 0.5*datavols->radial_boundary[0];
    for (int rr=1; rr<constants->nradii; rr++) {
        mean_radii[rr] = 0.5*(datavols->radial_boundary[rr]+datavols->radial_boundary[rr-1]);
    }
    checkCudaErrors( cudaMemcpyToSymbol( KERN_RADII, mean_radii, constants->nradii*sizeof(float) ) );
    delete[] mean_radii;
    /* checkCudaErrors( cudaMemcpyToSymbol( KERN_RADII, datavols->radial_boundary, constants->nradii*sizeof(float) ) ); */
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////// Bind Density to 3D Texture Array //////////////////////////////////////////////
    cudaExtent volExtent = make_cudaExtent(
            constants->size.x,
            constants->size.y,
            constants->size.z
            );
    memmgr.cudaMalloc3DArray(&d_dens_array[deviceid], &floatChannelDesc, volExtent);
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
    { // local scope
        // TODO: This depends strongly on the number of active beamlets for each beam and thus the
        // TODO Currently set statically in config.json
        // unsigned int actualBEVSize = make_uint3(max_diag_len + 2*(constants->penumbra) + buffer);
        // printf("dynamic max_bev size: X = %i,  Y = %i, Z = %i \n", constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
    }
    /* constants->max_rev_size = uint3{800,800,800}; */
    printf("  Full Convolution Memory Dimensions: %d x %d x %d\n",
            constants->max_rev_size.x,constants->max_rev_size.y,constants->max_rev_size.z);

    // Following memory operations are allocated on their assigned GPU

    // allocate intermediate data volume according to BEV extents // memsets occur later, just before use of each in kernel looping
    { // local scope
        int revSize = constants->max_rev_size.x*constants->max_rev_size.y*constants->max_rev_size.z*sizeof(float);
        memmgr.cudaMalloc( (void**)&device_data[deviceid].bevDens,      revSize );
        memmgr.cudaMalloc( (void**)&device_data[deviceid].bevTerma,     revSize );
    }

    // Cuda Array for terma texture object fetching
    cudaExtent calc_bbox_extent = make_cudaExtent(constants->calc_bbox_size.x, constants->calc_bbox_size.y, constants->calc_bbox_size.z);
    memmgr.cudaMalloc3DArray(&device_data[deviceid].term_Array, &floatChannelDesc, calc_bbox_extent);

    // allocate 3D CUDA Array to find with surface and texture objects
    { // local scope
        cudaExtent doseExtent = make_cudaExtent( constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
        memmgr.cudaMalloc3DArray(&(device_data[deviceid].dose_Array), &floatChannelDesc, doseExtent, cudaArraySurfaceLoadStore);
    }

    // generate texture/surface object for reading/writing terma and dose data in kernels
    makeTexObject<cudaArray>(&(device_data[deviceid].texTerma), device_data[deviceid].term_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);
    makeTexObject<cudaArray>(&(device_data[deviceid].texDose),  device_data[deviceid].dose_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);
    makeSurfObject(&(device_data[deviceid].surfDose), device_data[deviceid].dose_Array);
    //////////////////////////////////////////////////////////////////////////////////////////////

    if (debug_write) {
        memmgr.cudaMalloc( (void**)&device_dose_accum[deviceid], datavols->size_data * sizeof(float) );
        checkCudaErrors( cudaMemset( device_dose_accum[deviceid], 0, datavols->size_data * sizeof(float) ) );
    }

    if (verbose) {
        cudaMemInfo meminfo = query_device_memory();
        meminfo.gpuid = gpuid;
            meminfo.print_available();
    }
    std::cout << "Memory used on device "<<gpuid<<": "<<memmgr.memoryUsed()/GB<< " GB | available: "<<memmgr.memoryRemaining()/GB<<" GB"<<std::endl;
    if (timing) {
        timer.stop_print_time_elapsed("GPU mem allocation");
    }

    return 1;
}

PILLAR_GRID configure_beamlet_grid(
    BEAM const *beam,
    float3 rev_voxelsize,
    int verbose=0
) {
    return PILLAR_GRID();
}

// performs the convolution for each beam
int radconvolveTexture(
        MONO_KERNELS     *mono,
        CONSTANTS        *constants,
        BEAM             *device_beam_arr,
        int              nbeams,
        int              nrays,
        int              deviceid,      // non-CUDA device index (always start from 0, increment by 1)
        int              gpuid,         // CUDA device index (could skip ints for incompatible devices)
        int              num_unpack_streams,
        bool             verbose=false,
        bool             timing=false,
        bool             debugwrite=false
) {
    NVTX_PUSH_RANGE("2-Calc_Prepare", 2);
    cudaDeviceProp devProp;
    checkCudaErrors(cudaGetDeviceProperties(&devProp, gpuid));

    MemoryManager& memmgr = memory_managers[deviceid];

    float3 rev_voxelsize = {constants->rev_latspacing, constants->rev_longspacing, constants->rev_latspacing};

    if (TERMA_ONLY) {
        nrays = 1;
        /* constants->kernel_extent = 0.1; */
    }


    // set up arrays for dynamic GPU resource allocation
    // dependent on dimensions of BEV data per convolution ray
    dim3 tileBlock(TILE_DIM_X, TILE_DIM_Y, 1);
    dim3 rayGrid[nrays];
    dim3 conBlock[nrays];
    dim3 conGrid[nrays];
    unsigned int memsize[nrays];

    // Calculate cuda execution block/grid sizes
    unsigned int dataSize = constants->nvoxels();
    unsigned int calcDataSize = constants->bbox_nvoxels();
    unsigned int revSize = constants->max_rev_size.x * constants->max_rev_size.y * constants->max_rev_size.z;

    // create thread stream(s)
    cudaStream_t beamstream;
    cudaStreamCreate(&beamstream);
    {
        char stream_name[50];
        sprintf(stream_name, "BeamStream%d", deviceid);
        NVTX_NAME_STREAM(beamstream, stream_name);
    }

    // TODO: add as config file parameter
    float* h_unpacked_dose_arr[num_unpack_streams];
    float* d_unpacked_dose_arr[num_unpack_streams];
    cudaStream_t unpackstream_arr[num_unpack_streams];
    for (int n=0; n<num_unpack_streams; n++) {
        std::unique_ptr<float[]> temp(new float[calcDataSize*sizeof(float)]);
        sparsifymanager.push_memblock(temp); // add once-allocated memory blocks for use during sparsify/write-to-file
        checkCudaErrors( cudaHostAlloc((void**)&h_unpacked_dose_arr[n], calcDataSize*sizeof(float), 0) );
        memmgr.cudaMalloc((void**)&d_unpacked_dose_arr[n], calcDataSize*sizeof(float));
        cudaStreamCreate(&(unpackstream_arr[n]));
        char stream_name[50];
        sprintf(stream_name, "UnpackStream%d", n);
        NVTX_NAME_STREAM(unpackstream_arr[n], stream_name);
    }

    CudaTimer timer(beamstream);
    if (timing) { timer.start(); }

    /////////////// Declare Reusable Memory ///////////////////
    // allocate terma and dose data volumes // memsets performed at start of each stream/beam computation
    /* memmgr.cudaMalloc( (void**) &device_data[deviceid].terma, calcDataSize * sizeof(float) ); */

    float *d_fluence_map; // device storage of beamlet fluence weights/intensities
    memmgr.cudaMalloc( (void**) &d_fluence_map, device_beam_arr[0].fmap_size.x*device_beam_arr[0].fmap_size.y*sizeof(float) );

    float *h_beam_dose;  // host (pinned) memory for writing single beam dose volume to file
    checkCudaErrors( cudaHostAlloc( (void **) &h_beam_dose, dataSize * sizeof(float), 0) );
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////
    // Calculate memory requirements for all beams (vec[nbeams][nrays])
    /*---
    std::vector< std::vector<BLT_CONV_DATA> > conv_data(nbeams, std::vector<BLT_CONV_DATA>(nrays) );
    for (int bb=0; bb<nbeams; bb++) {


        for (int rr=0; rr<nrays; rr++) {
            findREV(&conv_data[bb][rr].rev,
                    constants,
                    &device_beam_arr[bb],
                    make_float3(0.f),
                    bPG.gridSize,
                    constants->get_theta_from_index(rr),
                    constants->get_phi_from_index(rr),
                    (verbose && extra_verbose)
                   );

        }
    }
    ---*/
    ///////////////////////////////////////////////////////////


    if (timing) { timer.reset_print_time_elapsed("Convolution setup"); }

    NVTX_POP_RANGE; // 2-Calc_Prepare
    NVTX_PUSH_RANGE("2-Calc_Beamloop", 2);

    uint3 max_actual_rev_size = {0};


    // for each beam assigned to this GPU stream
    for (int dc=0; dc<nbeams; dc++) {
        NVTX_PUSH_RANGE("3-Beam_Prepare", 3);


        BEAM* this_beam = &device_beam_arr[dc];
        std::cout<<"Starting: | device#"<<gpuid<<" beam#"<<dc<<" deviceid#"<<deviceid<<" | beam="<<*this_beam<<std::endl;
        if (timing) { timer.start(); }

        // collect per-beam memsets to increase stream concurrency
        // TODO: some memsets may be unnecessary if memcpy fills array with new vals anyways
        /* checkCudaErrors( cudaMemsetAsync(device_data[deviceid].terma,       0, calcDataSize * sizeof(float), beamstream) ); */
        // copy binary fluence mask to GPU
        checkCudaErrors( cudaMemcpyAsync(d_fluence_map, this_beam->fluence_map.data(), this_beam->fmap_size.x*this_beam->fmap_size.y*sizeof(float), cudaMemcpyHostToDevice, beamstream) );

        //////////////////////////////////////////////
        //////////////BEAMLET SPECIFIC////////////////
        // initialize beam-specific context grid
        PILLAR_GRID hPG{};

        // 1. get the number of beamlets.
        std::vector<int> temp_idx;
        for (int idx = 0; idx < this_beam->fmap_size.x * this_beam->fmap_size.y; idx++){
            if ( this_beam->fluence_map[ idx ] > 0 ){
                ++hPG.numBeamlets;
                temp_idx.push_back(idx);
            }
        }

        HEADER_BEAM beam_header;
        beam_header.beam_uid = this_beam->uid;
        beam_header.beam_specs = *this_beam;
        beam_header.N_beamlets = hPG.numBeamlets;

        int max_batch_size = hPG.numBeamlets; // max # beamlet contexts in any batch

        // prepare memory/vars for use later
        cudaArray* PackedBEVdose_Array;
        cudaTextureObject_t texPackedBEVDose;   // texture object created/destroyed for each beam - array mem is once initialized for all beams
        // TODO: try using float2 and float linear textures for getting pillarStartCoords and beamletAngles more efficiently in kernels (texture caching)
        int*    dpg_beamletIdx;
        float3* dpg_pillarStartCoords;
        float2* dpg_beamletAngles;
        float3* dpg_beamletIsocenters;
        checkCudaErrors( cudaHostAlloc( (void**) &hPG.beamletIdx, hPG.numBeamlets * sizeof( int ), 0 ) );
        checkCudaErrors( cudaHostAlloc( (void**) &hPG.pillarStartCoords, hPG.numBeamlets * sizeof( float3 ), 0 ) );
        checkCudaErrors( cudaHostAlloc( (void**) &hPG.beamletAngles, hPG.numBeamlets * sizeof( float2 ), 0 ) );
        checkCudaErrors( cudaHostAlloc( (void**) &hPG.beamletIsocenters, hPG.numBeamlets * sizeof( float3 ), 0 ) );
        memmgr.cudaMalloc( ( void**)& dpg_pillarStartCoords, max_batch_size * sizeof(float3) );
        memmgr.cudaMalloc( ( void**)& dpg_beamletAngles, max_batch_size * sizeof(float2) );
        memmgr.cudaMalloc( ( void**)& dpg_beamletIsocenters, max_batch_size * sizeof(float3) );
        memmgr.cudaMalloc( ( void**)& dpg_beamletIdx, max_batch_size * sizeof(int) );

        for (int ii=0; ii<temp_idx.size(); ii++) { hPG.beamletIdx[ii] = temp_idx[ii]; }

        // 2. get the central axis limits (anchors) of beamlets
        std::vector<float3> beamlet_start(hPG.numBeamlets);
        std::vector<float3> beamlet_end(hPG.numBeamlets);
        std::vector<float> beamlet_length(hPG.numBeamlets);
        float max_beamlet_length = 0.f;
        hPG.max_beamlet_size = float2{0.f, 0.f};
        for (int ii=0; ii<hPG.numBeamlets; ii++) {
            // get anchor positions for each beamlet
            calcBeamletAnchors(
                    beamlet_start[ii], beamlet_end[ii], hPG.beamletAngles[ii], hPG.beamletIsocenters[ii],
                    this_beam->source, this_beam->isocenter,
                    static_cast<unsigned int>(hPG.beamletIdx[ii]),
                    this_beam->beamlet_size, this_beam->fmap_size,
                    constants->voxel,
                    constants->start,
                    constants->calc_bbox_start,
                    constants->calc_bbox_size,
                    this_beam->azimuth, this_beam->zenith, this_beam->coll,
                    (int)verbose + (int)extra_verbose
                    );

            beamlet_length[ii] = length(beamlet_end[ii]-beamlet_start[ii]);
            if (beamlet_length[ii] > max_beamlet_length) { max_beamlet_length = beamlet_length[ii]; }

            float2 beamlet_diverge_size = this_beam->beamlet_size * length(beamlet_end[ii]-this_beam->source)/length(hPG.beamletIsocenters[ii]-this_beam->source);
            hPG.max_beamlet_size.x = max(hPG.max_beamlet_size.x, beamlet_diverge_size.x);
            hPG.max_beamlet_size.y = max(hPG.max_beamlet_size.y, beamlet_diverge_size.y);
        }

        // compute pillar size and ensure that it is an integer multiple of the rev_voxelsize
        float psize_long   = max_beamlet_length + 2.f*constants->kernel_extent + hPG.wallThickness*rev_voxelsize.y;
        // transverse size should use beamlet_size at largest end of beamlet (divergence effects)
        float2 psize_trans = hPG.max_beamlet_size + 2.f*constants->kernel_extent + hPG.wallThickness*rev_voxelsize.x;
        float3 expand = rev_voxelsize - make_float3(
                fmodf(psize_trans.x, rev_voxelsize.x),
                fmodf(psize_long,    rev_voxelsize.y),
                fmodf(psize_trans.y, rev_voxelsize.z) );
        // Also add buffer to each buffer so edges of pillar are not used (since these contain garbage from interpolation
        // during rotations)
        hPG.pillarSize = make_float3(psize_trans.x, psize_long, psize_trans.y) + expand + 2.f*float(hPG.pillarBuffer)*rev_voxelsize;
        hPG.pillarDims = make_int3(hPG.pillarSize/rev_voxelsize);

        // compute pillar limits to use in geometric transformations
        for (int ii=0; ii<hPG.numBeamlets; ii++) {
            float3 g_offset = make_float3(
                    -0.5f * (hPG.pillarSize.x + hPG.wallThickness*rev_voxelsize.x),
                    -0.5f * (hPG.pillarSize.y - beamlet_length[ii] + hPG.wallThickness*rev_voxelsize.y),
                    -0.5f * (hPG.pillarSize.z + hPG.wallThickness*rev_voxelsize.z) );
            // TODO fix rotation pivot to src here (may be unnecessary)
            g_offset = inverseRotateBeamAtOriginRHS(
                           inverseRotateBeamletAtOriginRHS(g_offset, hPG.beamletAngles[ii].x, hPG.beamletAngles[ii].y),
                       this_beam->azimuth, this_beam->zenith, this_beam->coll);
            hPG.pillarStartCoords[ii] = beamlet_start[ii] + g_offset;
            if (verbose && extra_verbose) {
                float3 pillarEndCoords = hPG.pillarStartCoords[ii] + hPG.pillarSize;
                printf(" ###%4i-beamlet %4i | start:(%8.3f, %8.3f, %8.3f); end:(%8.3f, %8.3f, %8.3f); length:%8.3f; diverge|azi:%7.2f, zen=%7.2f; pstart:(%8.3f, %8.3f, %8.3f); pend:(%8.3f, %8.3f, %8.3f)|\n",
                        ii, hPG.beamletIdx[ii],
                        beamlet_start[ii].x, beamlet_start[ii].y, beamlet_start[ii].z,
                        beamlet_end[ii].x, beamlet_end[ii].y, beamlet_end[ii].z,
                        beamlet_length[ii],
                        hPG.beamletAngles[ii].x * 180.f/PI,
                        hPG.beamletAngles[ii].y * 180.f/PI,
                        hPG.pillarStartCoords[ii].x, hPG.pillarStartCoords[ii].y, hPG.pillarStartCoords[ii].z,
                        pillarEndCoords.x, pillarEndCoords.y, pillarEndCoords.z
                      );
            }
        }

        if (verbose) {
            printf("beam %d.%d has %d beamlets\n", deviceid, dc, hPG.numBeamlets);
        }
        if (hPG.numBeamlets <= 0) {
            std::cout << set_color(COLOR::YELLOW) << "WARNING: beam " << deviceid<<"."<<dc<<" has 0 beamlets. Skipping" << set_color() << std::endl;
            // write empty beam entry to h5 file
            SRWorker::SRData data {};
            data.beam_header = beam_header;
            sparsifymanager.push(data);
            continue;
        }

        NVTX_POP_RANGE; // 3-Beam_Prepare
        NVTX_PUSH_RANGE("3-Beam_Batchloop", 3);

        int nbatches = 1;
        std::vector<int> batch_sizes = allocate_batch_sizes(hPG.numBeamlets, nbatches);
        int batch_memidx=0;
        for (int batchidx=0; batchidx<nbatches; batchidx++) {
            NVTX_PUSH_RANGE("4-Batch_Prepare", 4);

            // batch specific context grid
            PILLAR_GRID bPG{};
            //divide beamlets amongst batches
            bPG.beamletIdx = &hPG.beamletIdx[batch_memidx];
            bPG.pillarStartCoords = &hPG.pillarStartCoords[batch_memidx];
            bPG.beamletAngles = &hPG.beamletAngles[batch_memidx];
            bPG.beamletIsocenters = &hPG.beamletIsocenters[batch_memidx];
            bPG.numBeamlets = batch_sizes[batchidx];
            bPG.wallThickness = hPG.wallThickness;
            bPG.pillarBuffer = hPG.pillarBuffer;
            bPG.max_beamlet_size = hPG.max_beamlet_size;
            bPG.pillarDims = hPG.pillarDims;
            bPG.pillarSize = hPG.pillarSize;

            if (verbose && nbatches>1) {
                std::cout << "running beamlet batch "<<batchidx+1<< " of " << nbatches << " containing beamlets: "<< batch_memidx+1<<"->"<<batch_memidx+bPG.numBeamlets<<" of "<<hPG.numBeamlets<<std::endl;
            }
            batch_memidx += bPG.numBeamlets;
            if (batch_memidx > hPG.numBeamlets || max_batch_size < bPG.numBeamlets) {
                throw std::runtime_error("Problem with batching division. contact a developer.");
            }

            // 3. build the packed BEV allocation map
            bPG.numPillars.x = static_cast<int>( ceil(sqrt(float(bPG.numBeamlets))) );
            bPG.numPillars.y = static_cast<int>( ceil(float(bPG.numBeamlets)/bPG.numPillars.x) );

            // add another wall at the end of the array in each lateral dimension
            bPG.gridDims = make_int3(
                    bPG.numPillars.x*bPG.pillarDims.x + bPG.wallThickness,
                    bPG.pillarDims.y + bPG.wallThickness,
                    bPG.numPillars.y*bPG.pillarDims.z + bPG.wallThickness
                    );
            bPG.gridSize = make_float3(bPG.gridDims) * rev_voxelsize;

            if (verbose && extra_verbose) {
                printf("\n");
                printf("  # Active Beamlets:     %d\n", bPG.numBeamlets);
                printf("  Max Beamlet Size:      %6.2f x %6.2f\n", bPG.max_beamlet_size.x, bPG.max_beamlet_size.y);
                printf("  Iso Beamlet Size:      %6.2f x %6.2f\n", this_beam->beamlet_size.x, this_beam->beamlet_size.y);
                printf("  # Pillars:             %d x %d\n", bPG.numPillars.x, bPG.numPillars.y);
                printf("  Pillar Dims:           %d x %d x %d\n", bPG.pillarDims.x, bPG.pillarDims.y, bPG.pillarDims.z);
                printf("  Pillar Size [cm]:      (%6.2f, %6.2f, %6.2f)\n", bPG.pillarSize.x, bPG.pillarSize.y, bPG.pillarSize.z);
                printf("  PillarGrid Dims:       %d x %d x %d\n", bPG.gridDims.x, bPG.gridDims.y, bPG.gridDims.z);
                printf("  PillarGrid Size [cm]:  (%6.2f, %6.2f, %6.2f)\n\n", bPG.gridSize.x, bPG.gridSize.y, bPG.gridSize.z);
            }
            //////////////END BEAMLET SPECIFIC////////////////
            //////////////////////////////////////////////////

            dim3 rayBlock;
            constants->ss_factor = 1; // XXX: Need to fix to include support for >1
            if (constants->ss_factor>1) {
                uint nsubvoxels = powf(constants->ss_factor, 3);
                uint rayThreadLimit = 256;
                uint tilesize = floorf(sqrtf(rayThreadLimit/nsubvoxels));
                rayBlock = dim3{nsubvoxels, tilesize, tilesize};
            } else {
                rayBlock = dim3{1, 32, 3};
            }
            int raySharedMem = rayBlock.x*rayBlock.y*rayBlock.z*sizeof(float);

            // holds dimensions and data for the active volume when viewed at each convolution angle
            // used for dynamic GPU memory and resource allocation
            REV_DATA rev[nrays];

            // find the dimensions of the transformed data in BEV/NVB coordinate system for all rays
            // rev->min_coords and rev->max_coords will be an expansion of pbev_lim_min and pbev_lim_max based on REV rotation of pillar_grid
            try {
                for (int rr=0; rr<nrays; rr++) {
                    float theta = constants->get_theta_from_index(rr);
                    float phi = constants->get_phi_from_index(rr);
                    if (TERMA_ONLY) {
                        theta = 0.f;
                        phi = 0.f;
                    }
                    // compute angles from idx
                    if (extra_verbose) { printf(" Ray #%d\n", rr+1); }
                    findREV(&rev[rr],
                            constants,
                            this_beam,
                            make_float3(0.f),
                            bPG.gridSize,
                            theta,
                            phi,
                            (verbose && extra_verbose)
                           );

                    // PAY ATTENTION: To keep memory accesses coalesced where possible coordinates and data layout axes are not the same
                    //     as such, kernel launch param axes also do not match. coordinate axes XYZ are mapped to device memory layout axes ZXY
                    rayGrid[rr] = make_uint3(
                            rev[rr].size.z,
                            static_cast<unsigned int>(ceilf(static_cast<float>(rev[rr].size.x)/rayBlock.y)),
                            static_cast<unsigned int>(ceilf(static_cast<float>(rev[rr].size.y)/rayBlock.z))
                            );

                    // x    # samples along conv. ray direction
                    // y|z: # parallel rays in REV
                    // TODO: be careful when rev[rr].size.x approaches #threads-per-block limit (1024 for CC3.0+)
                    //       In most cases this should be below 300-400 though
                    if (rev[rr].size.x > 1024) { throw std::runtime_error("size of rev (x-dim) exceeds maximum #threads allowed. Contact dev to resolve."); }
                    int overwrite = 2; // buffer beyond used rev block, where convolve kernel will also write zeros before proceeding
                    conBlock[rr].x = rev[rr].size.x + overwrite;
                    conGrid[rr].y  = rev[rr].size.y + overwrite;
                    conGrid[rr].z  = rev[rr].size.z + overwrite;
                    conBlock[rr].y = conBlock[rr].z = conGrid[rr].x = 1;
                    // dynamic alloc. shared mem size
                    memsize[rr] = 2 * (conBlock[rr].x - overwrite) * sizeof(float);

                    if (verbose && extra_verbose) {
                        printf("Terma Supersampling: %s\n", (constants->ss_factor>1)?"enabled":"disabled");
                        printf(" ### rayGrid:  %3d x %3d x %3d  |  rayBlock:  %3d x %3d x %3d | Shared: %d bytes\n",
                                rayGrid[rr].x, rayGrid[rr].y, rayGrid[rr].z,
                                rayBlock.x, rayBlock.y, rayBlock.z,
                                raySharedMem
                              );
                        printf(" conGrid:   %3d x %3d x %3d  |  conBlock:   %3d x %3d x %3d\n\n",
                                conGrid[rr].x, conGrid[rr].y, conGrid[rr].z,
                                conBlock[rr].x, conBlock[rr].y, conBlock[rr].z);
                    }

                    max_actual_rev_size.x = max(max_actual_rev_size.x, rev[rr].size.x);
                    max_actual_rev_size.y = max(max_actual_rev_size.y, rev[rr].size.y);
                    max_actual_rev_size.z = max(max_actual_rev_size.z, rev[rr].size.z);
                }
            } catch (const ArraySizeError& e) {
                if (nbatches > 100) {
                    throw std::runtime_error("Not enough GPU memory for the selected quality settings. Please reduce quality and retry.");
                }
                std::cout << set_color(COLOR::YELLOW)<<"Convolution array size for ray (az:"<<e.theta*180.f/PI<<"|zen:"<<e.phi*180.f/PI<<" deg) is larger than allocated memory"
                    << " ("<<e.oldsize.x<<","<<e.oldsize.y<<","<<e.oldsize.z<<") > ("<<e.rev->size.x<<","<<e.rev->size.y<<","<<e.rev->size.z<<"). "
                    << "Incrementing batch count by 1" << set_color() << std::endl;
                batch_memidx = 0;
                batchidx = -1;
                nbatches++;
                batch_sizes = allocate_batch_sizes(hPG.numBeamlets, nbatches);
                continue;
                /* if (!debugwrite) { throw e; } */
            }


            // Prepare kernel launch parameters
            dim3 packedGrid = dim3(
                    static_cast<unsigned int>(ceilf(static_cast<float>(bPG.gridDims.y)/tileBlock.x)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(bPG.gridDims.z)/tileBlock.y)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(bPG.gridDims.x)/tileBlock.z))
                    );
            if (verbose && extra_verbose) {
                printf(" ### packedGrid:   %3d x %3d x %3d  |  packedBlock:   %3d x %3d x %3d\n",
                        packedGrid.x, packedGrid.y, packedGrid.z,
                        tileBlock.x, tileBlock.y, tileBlock.z
                      );
            }


            // Move pillar_grid data to device
            checkCudaErrors( cudaMemcpyAsync( dpg_beamletIdx, bPG.beamletIdx, bPG.numBeamlets * sizeof(int), cudaMemcpyHostToDevice, beamstream ));
            checkCudaErrors( cudaMemcpyAsync( dpg_pillarStartCoords, bPG.pillarStartCoords, bPG.numBeamlets * sizeof(float3), cudaMemcpyHostToDevice, beamstream ));
            checkCudaErrors( cudaMemcpyAsync( dpg_beamletAngles, bPG.beamletAngles, bPG.numBeamlets * sizeof(float2), cudaMemcpyHostToDevice, beamstream ));
            checkCudaErrors( cudaMemcpyAsync( dpg_beamletIsocenters, bPG.beamletIsocenters, bPG.numBeamlets * sizeof(float3), cudaMemcpyHostToDevice, beamstream ));

            // TODO we can probably deterimine the max size needed here when we choose max_rev_size based on number of beamlets, max beamlet width and kernel_extent
            memmgr.cudaMalloc( (void**) &device_data[deviceid].dose, bPG.pillar_grid_nvoxels() * sizeof(float) );
            checkCudaErrors( cudaMemsetAsync(device_data[deviceid].dose, 0, bPG.pillar_grid_nvoxels() * sizeof(float), beamstream) );

            if (timing) { timer.restart_print_time_elapsed("Beam setup"); }

            // ---------------------------------------------------------------

            NVTX_POP_RANGE;
            NVTX_PUSH_RANGE("4-Batch_Convloop", 4);

            int debug_ray = constants->ntheta/2+1;
            for (int ray_idx=0; ray_idx<nrays; ray_idx++) {
                // get kernel rotation angles once instead of for each kernel thread launch
                int kern_wt_idx = constants->get_kernel_theta_index(ray_idx);
                float kern_theta = constants->get_theta_from_index(ray_idx);
                float kern_phi = constants->get_phi_from_index(ray_idx);
                if (TERMA_ONLY) {
                    kern_theta = 0.f;
                    kern_phi = 0.f;
                }

                /*---
                dim3 siddonBlock;
                dim3 siddonGrid;
                if (constants->ss_factor>1) {
                    uint nsubvoxels = powf(constants->ss_factor, 3);
                    uint tilesize = floorf(sqrtf(256/nsubvoxels));
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
                if (verbose && extra_verbose) {
                    printf("Terma Supersampling: %s\n", (constants->ss_factor>1)?"enabled":"disabled");
                    printf(" ### siddonGrid:  %3d x %3d x %3d  |  siddonBlock:  %3d x %3d x %3d | Shared: %d bytes\n",
                            siddonGrid.x, siddonGrid.y, siddonGrid.z,
                            siddonBlock.x, siddonBlock.y, siddonBlock.z,
                            siddonSharedMem
                          );
                }

                float* accum_terma;
                if (debugwrite && ray_idx==0) {
                    memmgr.cudaMalloc( (void**)&accum_terma, calcDataSize*sizeof(float) );
                }

                for (int bidx=0; bidx<bPG.numBeamlets; bidx++) {
                    int beamletnum = bPG.beamletIdx[bidx];
                    // printf("processing ray %d, beamlet %d\n", ray_idx, beamletnum);
                    checkCudaErrors( cudaMemsetAsync(device_data[deviceid].terma,       0, calcDataSize * sizeof(float), beamstream) );
                    cudaSiddon<<<siddonGrid, siddonBlock, siddonSharedMem, beamstream>>>(
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
                        this_beam->beamlet_size,
                        this_beam->fmap_size,
                        this_beam->azimuth,
                        this_beam->zenith,
                        this_beam->coll,
                        mono->nkernels,
                        texDens[deviceid],
                        texSpectrum[deviceid],
                        beamletnum
                        );

                    // write debug output (terma) to file
                    if (debugwrite && ray_idx==0 ) {
                        if ( extra_debug ) {
                            // save terma to file
                            checkCudaErrors( cudaStreamSynchronize(beamstream) );
                            float* out_array = new float[calcDataSize * sizeof(float) ];
                            memset((void*)out_array, 0, calcDataSize * sizeof(float));
                            checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].terma, calcDataSize*sizeof(float), cudaMemcpyDeviceToHost));

                            char checkchar[50];
                            sprintf(checkchar,"terma-d%d-B%d-b%d", gpuid, dc, beamletnum);
                            write_debug_data<float>( out_array, constants->calc_bbox_size, checkchar, verbose);
                            delete [] out_array;
                        }

                        // accumulate terma
                        cudaSum<<<(int)ceilf(calcDataSize/1024.f), 1024>>>(accum_terma, device_data[deviceid].terma, calcDataSize);
                    }

                    // copy terma volume to 3D cuda array, and bind to texture with local scope
                    // // term_array[_threadid] allocated outside of beam loop
                    { // local scope
                        cudaExtent calc_bbox_extent = make_cudaExtent(constants->calc_bbox_size.x, constants->calc_bbox_size.y, constants->calc_bbox_size.z);
                        cudaMemcpy3DParms CopyParams = {0};
                        CopyParams.srcPtr   = make_cudaPitchedPtr((void*)device_data[deviceid].terma, calc_bbox_extent.width*sizeof(float), calc_bbox_extent.width, calc_bbox_extent.height);
                        CopyParams.dstArray = device_data[deviceid].term_Array;
                        CopyParams.extent   = calc_bbox_extent;
                        CopyParams.kind     = cudaMemcpyDeviceToDevice;
                        checkCudaErrors( cudaMemcpy3DAsync(&CopyParams, beamstream) );
                    }
                }
                ---*/
                {
                    cudaBeamletRaytrace <<< rayGrid[ray_idx], rayBlock, raySharedMem, beamstream >>> (
                            device_data[deviceid].bevDens,      // BEV output
                            device_data[deviceid].bevTerma,     // BEV output
                            this_beam->source,
                            this_beam->beamlet_size,
                            this_beam->azimuth,
                            this_beam->zenith,
                            this_beam->coll,
                            make_float3(bPG.gridDims),
                            make_float3(bPG.pillarDims),
                            bPG.numBeamlets,
                            bPG.wallThickness,
                            bPG.numPillars,
                            dpg_beamletIdx,
                            dpg_beamletAngles,
                            dpg_pillarStartCoords,
                            dpg_beamletIsocenters,
                            kern_theta, kern_phi,                // ray direction
                            rev[ray_idx].min_coords,             // bev box start indices // but i strongly doubt this is GCS coordinates
                            rev[ray_idx].size,                   // BEV volume dims
                            constants->max_rev_size,
                            constants->start,
                            constants->voxel,
                            rev_voxelsize,
                            make_float3(constants->calc_bbox_start),
                            texDens[deviceid],                   // XYZ input
                            device_data[deviceid].texTerma,     // XYZ input

                            // raytracing/terma args
                            d_fluence_map,
                            make_float3(constants->size),
                            constants->calc_bbox_size,
                            constants->beamhard_correct,
                            this_beam->direction,
                            this_beam->isocenter,
                            this_beam->sad,
                            this_beam->fmap_size,
                            mono->nkernels,
                            texSpectrum[deviceid]
                                );
                    if (debugwrite && extra_debug && ray_idx <= debug_ray) {
                        checkCudaErrors( cudaStreamSynchronize(beamstream) );
                        float* out_array = new float[revSize * sizeof(float) ];
                        memset((void*)out_array, 0, revSize * sizeof(float));
                        checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].bevDens, revSize*sizeof(float), cudaMemcpyDeviceToHost));

                        char checkchar[40];
                        sprintf(checkchar,"debug_dens-d%d-B%d-b%d-r%d", gpuid, dc, batchidx, ray_idx);
                        write_debug_data<float>(out_array, constants->max_rev_size, checkchar, true );
                        checkCudaErrors( cudaStreamSynchronize(beamstream) );
                        /////////////////////////////////////////////////////
                        memset((void*)out_array, 0, revSize * sizeof(float));
                        checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].bevTerma, revSize*sizeof(float), cudaMemcpyDeviceToHost));

                        sprintf(checkchar,"debug_terma-d%d-B%d-b%d-r%d", gpuid, dc, batchidx, ray_idx);
                        write_debug_data<float>(out_array, constants->max_rev_size, checkchar, true );
                        checkCudaErrors( cudaStreamSynchronize(beamstream) );

                        delete [] out_array;
                        checkCudaErrors(cudaStreamSynchronize(beamstream));
                    }
                }

                /*---
                if (debugwrite && ray_idx==0) {
                    cudaStreamSynchronize(beamstream);
                    float *out_array = new float[calcDataSize];
                    cudaMemcpy(out_array, accum_terma, calcDataSize*sizeof(float), cudaMemcpyDeviceToHost);
                    char checkchar[50];
                    sprintf(checkchar,"accum_terma-d%d-B%d", gpuid, dc);
                    write_debug_data<float>(out_array, constants->calc_bbox_size, checkchar, verbose);
                    delete[] out_array;
                    cudaFree(accum_terma);
                }
                ---*/

                // perform dose calcluation (CCCS) w/ heterogeneity correction in REV volume
                PackRowConvolve <<< conGrid[ray_idx], conBlock[ray_idx], memsize[ray_idx], beamstream >>> (
                        device_data[deviceid].bevDens,      // BEV-trans input
                        device_data[deviceid].bevTerma,     // BEV-trans input
                        device_data[deviceid].surfDose,     // BEV dose output
                        (float)kern_wt_idx,                  // kernel theta index (for sampling texKern)
                        rev[ray_idx].size,                   // BEV volume dims
                        constants->max_rev_size,
                        constants->rev_longspacing,
                        constants->nradii,
                        constants->ntheta,
                        constants->nphi,
                        texKern[deviceid],                    // kernel weights
                        TERMA_ONLY
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
                    sprintf(checkchar,"debug_dose-d%d-B%d-b%d-r%d", gpuid, dc, batchidx, ray_idx);
                    write_debug_data<float>(out_array, constants->max_rev_size, checkchar, true);

                    delete [] out_array;
                    checkCudaErrors(cudaStreamSynchronize(beamstream));
                }

                /* transform packed REV dose coefficients from the previous convolution back to BEV system then */
                /* perform element-by-element sum, accumulating over all convolution directions */
                PackedREVtoBEVdose <<< packedGrid, tileBlock, 0, beamstream >>> (
                        device_data[deviceid].dose,    // beamlet-packed dose array in BEV orientation
                        device_data[deviceid].texDose, // packed dose array embedded in REV bounding box
                        kern_theta, kern_phi,          // convolution direction
                        rev[ray_idx].min_coords,                    // REV volume limit coords in XYZ coord system
                        rev_voxelsize,
                        bPG.gridDims
                        );

                if (debugwrite && extra_debug && ray_idx <= debug_ray) {
                    int pgSize = bPG.pillar_grid_nvoxels();
                    checkCudaErrors( cudaStreamSynchronize(beamstream) );
                    float* out_array = new float[pgSize * sizeof(float) ];
                    memset((void*)out_array, 0, pgSize * sizeof(float));
                    checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].dose, pgSize*sizeof(float), cudaMemcpyDeviceToHost));

                    char checkchar[40];
                    sprintf(checkchar,"debug_bev_packeddose-d%d-B%d-b%d-r%d", gpuid, dc, batchidx, ray_idx);
                    uint3 pbev_write_dims = make_uint3(bPG.gridDims.y, bPG.gridDims.z, bPG.gridDims.x);
                    write_debug_data<float>(out_array, pbev_write_dims, checkchar, true);

                    delete [] out_array;
                }
            }
            cudaStreamSynchronize(beamstream);
            getLastCudaError("Kernel execution failed. Likely due to insufficient GPU memory. Reduce quality and retry.");

            // see Konvolve code + Papanikolaou and Mackie 1993

            if ( debugwrite) {
                int pgSize = bPG.pillar_grid_nvoxels();
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
                float* out_array = new float[pgSize * sizeof(float) ];
                memset((void*)out_array, 0, pgSize * sizeof(float));
                checkCudaErrors(cudaMemcpy((void*)out_array, (void*)device_data[deviceid].dose, pgSize*sizeof(float), cudaMemcpyDeviceToHost));

                char checkchar[40];
                sprintf(checkchar,"debug_bev_packeddose-d%d-B%d-b%d", gpuid, dc, batchidx);
                uint3 pbev_write_dims = make_uint3(bPG.gridDims.y, bPG.gridDims.z, bPG.gridDims.x);
                write_debug_data<float>(out_array, pbev_write_dims, checkchar, true);

                // also write density in packed bev
                cudaExtent doseExtent = make_cudaExtent( constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
                cudaMemcpy3DParms CopyParams = {0};
                CopyParams.srcPtr   = make_cudaPitchedPtr(device_data[deviceid].bevDens, sizeof(float)*doseExtent.width, doseExtent.width, doseExtent.height);
                CopyParams.dstArray = device_data[deviceid].dose_Array;
                CopyParams.extent   = doseExtent;
                CopyParams.kind     = cudaMemcpyDeviceToDevice;
                checkCudaErrors(cudaMemcpy3D(&CopyParams));

                float* d_out_array;
                memmgr.cudaMalloc((void**)&d_out_array, pgSize*sizeof(float));
                memset(out_array, 0, pgSize * sizeof(float));
                float kern_theta = constants->get_theta_from_index(nrays-1);
                float kern_phi = constants->get_phi_from_index(nrays-1);
                if (TERMA_ONLY) {
                    kern_theta = 0.f;
                    kern_phi = 0.f;
                }
                checkCudaErrors( cudaMemsetAsync(d_out_array, 0, pgSize*sizeof(float), beamstream) );
                checkCudaErrors( cudaStreamSynchronize(beamstream) );
                PackedREVtoBEVdose <<< packedGrid, tileBlock, 0, beamstream >>> (
                        d_out_array,    // beamlet-packed dose array in BEV orientation
                        device_data[deviceid].texDose, // packed dose array embedded in REV bounding box
                        kern_theta, kern_phi,          // convolution direction
                        rev[nrays-1].min_coords,                    // REV volume limit coords in XYZ coord system
                        rev_voxelsize,
                        bPG.gridDims
                        );
                checkCudaErrors( cudaMemcpy(out_array, d_out_array, pgSize*sizeof(float), cudaMemcpyDeviceToHost) );
                sprintf(checkchar, "debug_bev_packeddens-d%d-B%d-b%d", gpuid, dc, batchidx);
                write_debug_data<float>(out_array, pbev_write_dims, checkchar, true);

                delete [] out_array;
                memmgr.cudaFree(d_out_array);
            }

            if (timing) {
                timer.stop();
                float cudatime = timer.time_elapsed();
                printf("-- %d.%d - Kernel Convolve Time: %4.3f msec (avg time: %4.3f msec per ray (of %d rays) --\n", gpuid, dc, cudatime, cudatime/nrays, nrays);
                printf("\n"); // blankline between beams of same device
            }


            ////////// Unpack Pillars from BEV storage, sparsify, and write to disk
            // copy output of PackedREVtoBEVdose to cudaArray and attach texture object
            //TODO: move this into init method
            cudaChannelFormatDesc floatChannelDesc = cudaCreateChannelDesc<float>();
            cudaExtent packedArrayExtent = make_cudaExtent(bPG.gridDims.y, bPG.gridDims.z, bPG.gridDims.x);
            memmgr.cudaMalloc3DArray(&PackedBEVdose_Array, &floatChannelDesc, packedArrayExtent);
            // copy to cudaArray
            cudaMemcpy3DParms copyParams = {0};
            copyParams.srcPtr = make_cudaPitchedPtr((void*)device_data[deviceid].dose, packedArrayExtent.width*sizeof(float), packedArrayExtent.width, packedArrayExtent.height);
            copyParams.dstArray = PackedBEVdose_Array;
            copyParams.extent = packedArrayExtent;
            copyParams.kind = cudaMemcpyDeviceToDevice;
            checkCudaErrors( cudaMemcpy3DAsync(&copyParams, beamstream) );
            // attach to texture object
            makeTexObject<cudaArray>(&texPackedBEVDose, PackedBEVdose_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);

            // determine kernel launch params
            dim3 unpackBlock = tileBlock;
            dim3 unpackGrid = dim3(
                    static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.x)/unpackBlock.x)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.y)/unpackBlock.y)),
                    static_cast<unsigned int>(ceilf(static_cast<float>(constants->calc_bbox_size.z)/unpackBlock.z))
                    );

            cudaStreamSynchronize(beamstream);

            NVTX_POP_RANGE; // 4-Batch_Convloop
            NVTX_PUSH_RANGE("4-Batch_Unpack", 4);

            Timer timer_sparsify;
            if (timing) { timer_sparsify.start(); }

            // array properties needed for sparsify
            ArrayProps props {};
            props.size = constants->size;
            props.crop_size = constants->calc_bbox_size;
            props.crop_start = constants->calc_bbox_start;

            Timer timer_unpack;
            Timer timer_wait;
            Timer timer_sr;
            Timer timer_sync;
            int runningPillarIndex = 0;
            while (runningPillarIndex < bPG.numBeamlets) {
                int batchsize = min(num_unpack_streams, bPG.numBeamlets-runningPillarIndex);
                for (int streamid=0; streamid<batchsize; streamid++) {
                    if (timing) { timer_unpack.start(); }
                    int pillarIndex = runningPillarIndex + streamid;

                    // enable concurrent cuda streams
                    float* h_unpacked_dose_pinned = h_unpacked_dose_arr[streamid];
                    float* d_unpacked_dose = d_unpacked_dose_arr[streamid];
                    cudaStream_t& unpackstream = unpackstream_arr[streamid];

                    // re-initialize d_unpacked_dose
                    memset(h_unpacked_dose_pinned, 0, calcDataSize*sizeof(float));
                    checkCudaErrors( cudaMemsetAsync(d_unpacked_dose, 0, calcDataSize*sizeof(float), unpackstream) );

                    // Unpack Pillar from texpackedBEVDose
                    UnpackBEVDosePillar <<< unpackGrid, unpackBlock, 0, unpackstream >>>(
                            d_unpacked_dose,
                            texPackedBEVDose,
                            this_beam->sad,
                            this_beam->source,
                            this_beam->azimuth,
                            this_beam->zenith,
                            this_beam->coll,
                            constants->start,
                            constants->voxel,
                            make_float3(constants->calc_bbox_start),
                            constants->calc_bbox_size,
                            rev_voxelsize,
                            pillarIndex,
                            pillarIndex % bPG.numPillars.x,
                            pillarIndex / bPG.numPillars.x,
                            bPG.pillarDims,
                            bPG.wallThickness,
                            bPG.pillarBuffer,
                            bPG.pillarStartCoords[pillarIndex],
                            bPG.beamletAngles[pillarIndex]
                            );
                    getLastCudaError("Kernel execution failed");

                    // accumulate dose in one volume to check for artifacts
                    if (debugwrite) {
                        // determine kernel launch params
                        dim3 accumBlock(devProp.maxThreadsPerBlock/4,1,1);
                        unsigned int xdim = static_cast<unsigned int>(ceilf(static_cast<float>(constants->bbox_nvoxels())/accumBlock.x));
                        dim3 accumGrid = dim3(xdim, 1, 1);
                        cudaSumSubArray <<< accumGrid, accumBlock, 0, unpackstream>>>(device_dose_accum[deviceid], d_unpacked_dose, constants->calc_bbox_start, constants->calc_bbox_size, constants->size);
                        getLastCudaError("Kernel execution failed");
                    }

                    // memcpy d_unpacked_dose to host
                    checkCudaErrors( cudaMemcpyAsync(h_unpacked_dose_pinned, d_unpacked_dose, calcDataSize*sizeof(float), cudaMemcpyDeviceToHost, unpackstream) );
                    if (timing) { timer_unpack.stop(); }
                } // end batch loop (1)

                /* start next loop */
                // TODO: consider processing first ready stream instead of first indexed stream
                /* std::vector<bool> pendingstreams(batchsize, true); */
                /* int streamid = -1; */
                for (int streamid=0; streamid<batchsize; ++streamid) {
                /* while (std::any_of(pendingstreams.cbegin(), pendingstreams.cend(), [](bool i){ return i; })) { */
                    /* streamid = (++streamid)%batchsize; */

                    NVTX_PUSH_RANGE("5.1-InlinePostProcessing", 5);
                    int pillarIndex = runningPillarIndex + streamid;
                    int beamletIndex = bPG.beamletIdx[pillarIndex];
                    cudaStream_t& unpackstream = unpackstream_arr[streamid];

                    // SYNC happens here until a memblock can be obtained
                    NVTX_PUSH_RANGE("5.2-WaitForSparsifyQueue", 6);
                    if (timing) { timer_wait.start(); }
                    std::unique_ptr<float[]> h_unpacked_dose = sparsifymanager.get_memblock();
                    /* if (!h_unpacked_dose) { */
                    /*     continue; */
                    /* } */
                    if (timing) { timer_wait.stop(); }
                    NVTX_POP_RANGE;

                    NVTX_PUSH_RANGE("5.2-WaitForContextUnpack", 7);
                    if (timing) { timer_sync.start(); }
                    checkCudaErrors( cudaStreamSynchronize(unpackstream) ); // block until this stream's data is copied to pinned host mem
                    if (timing) { timer_sync.stop(); }
                    NVTX_POP_RANGE;

                    NVTX_PUSH_RANGE("5.2-CopyDoseToSparsifyQueue", 8);
                    if (timing) { timer_sr.start(); }

                    float* h_unpacked_dose_pinned = h_unpacked_dose_arr[streamid];
                    memcpy(h_unpacked_dose.get(), h_unpacked_dose_pinned, calcDataSize*sizeof(float));

                    if (debugwrite && (beamletIndex == 44 || beamletIndex == 45 || beamletIndex == 54 || beamletIndex == 55)) {
                        char outname[100];
                        sprintf(outname, "unpacked_dose-d%d-B%d-ii%04d-b%04d", gpuid, dc, pillarIndex+batchidx*max_batch_size, beamletIndex );
                        write_debug_data<float>(h_unpacked_dose.get(), constants->calc_bbox_size, outname, false);
                    }

                    // collect header data
                    HEADER_BEAMLET beamlet_header;
                    beamlet_header.beamlet_uid = beamletIndex;

                    // add processing task to queue: sparsify/write to disk thread
                    SRWorker::SRData data {};
                    data.dense_array = std::move(h_unpacked_dose);
                    data.props = props;
                    data.beamlet_header = std::move(beamlet_header);
                    data.beam_header = beam_header;
                    data.roi_list = &roi_list;

                    sparsifymanager.push(data);
                    // host memory for h_unpacked_dose is free'd inside of worker thread after each job is finished
                    if (timing) { timer_sr.stop(); }
                    NVTX_POP_RANGE;
                    NVTX_POP_RANGE;

                    /* pendingstreams[streamid] = false; */
                } // end batch loop (2)
                runningPillarIndex += batchsize;
            } // end unpack loop

            if (timing) { timer_unpack.stop_print_time_elapsed("Unpacking Pillars"); }
            if (timing) { timer_sync.stop_print_time_elapsed("Syncing after unpacking dose"); }
            if (timing) { timer_wait.stop_print_time_elapsed("Waiting for available queue slot"); }
            if (timing) { timer_sr.stop_print_time_elapsed("Sending data to SRWorkers"); }
            if (timing) { timer_sparsify.stop_print_time_elapsed("Write Sparse Beamlet Coefficients"); }

            NVTX_POP_RANGE; // 4-Batch_Convloop
            NVTX_PUSH_RANGE("4-Batch_Cleanup", 4); // 4-Batch_Convloop

            // TODO: can we take this out of the loop?
            memmgr.cudaFree(device_data[deviceid].dose);

            // Move this out of the beam loop - memset to reinitialize instead
            checkCudaErrors( cudaDestroyTextureObject(texPackedBEVDose) );
            memmgr.cudaFreeArray(PackedBEVdose_Array);

            NVTX_POP_RANGE; // 4-Batch_Cleanup

        } // End of beamlet batching

        NVTX_POP_RANGE; // 3-Beam_Batchloop
        NVTX_PUSH_RANGE("3-Beam_Cleanup", 3);

        if (debugwrite){
            float* h_device_dose_accum = new float[constants->nvoxels()]();
            checkCudaErrors( cudaMemcpy(h_device_dose_accum, device_dose_accum[deviceid], constants->nvoxels()*sizeof(float), cudaMemcpyDeviceToHost) );
            char outname[100];
            sprintf(outname, "accumulated_dose-d%d-B%d", gpuid, dc);
            write_debug_data<float>(h_device_dose_accum, constants->size, outname, true);
            delete [] h_device_dose_accum;
        }



//////////////////////////////////////////////////////////////////////////////////////////

        // TODO: Move this outside of loop - reinit between beams
        memmgr.cudaFree(dpg_beamletIdx);
        memmgr.cudaFree(dpg_pillarStartCoords);
        memmgr.cudaFree(dpg_beamletAngles);
        memmgr.cudaFree(dpg_beamletIsocenters);
        checkCudaErrors( cudaFreeHost(hPG.beamletIdx) );
        checkCudaErrors( cudaFreeHost(hPG.pillarStartCoords) );
        checkCudaErrors( cudaFreeHost(hPG.beamletAngles) );
        checkCudaErrors( cudaFreeHost(hPG.beamletIsocenters) );

        NVTX_POP_RANGE; // 3-Beam_Cleanup
    }
    NVTX_POP_RANGE; // 2-Calc_Beamloop
    NVTX_PUSH_RANGE("2-Calc_Cleanup", 2);
    if (timing) { timer.start(); }

    // cleanup thread memory
    memmgr.cudaFree(d_fluence_map);
    checkCudaErrors( cudaFreeHost(h_beam_dose) );
    for (int n=0; n<num_unpack_streams; n++) {
        checkCudaErrors( cudaFreeHost(h_unpacked_dose_arr[n]) );
        memmgr.cudaFree(d_unpacked_dose_arr[n]);
    }

    if (timing) { timer.stop_print_time_elapsed("Convolution cleanup"); }

    if (debugwrite){
        float* h_device_dose_accum = new float[constants->nvoxels()]();
        checkCudaErrors( cudaMemcpy(h_device_dose_accum, device_dose_accum[deviceid], constants->nvoxels()*sizeof(float), cudaMemcpyDeviceToHost) );
        char outname[100];
        sprintf(outname, "accumulated_dose-d%d", gpuid);
        write_debug_data<float>(h_device_dose_accum, constants->size, outname, true);
        delete [] h_device_dose_accum;
    }

    if (verbose) {
        printf("Maximum REV size encountered was: (%d, %d, %d)\n", max_actual_rev_size.x, max_actual_rev_size.y, max_actual_rev_size.z);
        printf("Maximum GPU Memory usage (dev:0) was: %0.3f MB\n", memory_managers[0].peakMemoryUsed()/1024.f/1024.f);
    }

    NVTX_POP_RANGE; // 2-Calc_Cleanup
    return true;
}

// free GPU resources for each fork
int
freeCudaTexture(int ndevices, int *gpuid_arr, bool debugwrite)
{
    for (int deviceid=0; deviceid<ndevices; deviceid++) {
        MemoryManager& memmgr = memory_managers[deviceid];
        int gpuid = gpuid_arr[deviceid];
        checkCudaErrors( cudaSetDevice(gpuid) );

        memmgr.cudaFree(device_data[deviceid].bevDens);
        memmgr.cudaFree(device_data[deviceid].bevTerma);
        checkCudaErrors(cudaDestroySurfaceObject(device_data[deviceid].surfDose));
        checkCudaErrors(cudaDestroyTextureObject(device_data[deviceid].texDose));
        checkCudaErrors(cudaDestroyTextureObject(device_data[deviceid].texTerma));
        memmgr.cudaFreeArray(device_data[deviceid].dose_Array);
        memmgr.cudaFreeArray(device_data[deviceid].term_Array);
        /* memmgr.cudaFree(device_data[deviceid].terma); */

        checkCudaErrors(cudaDestroyTextureObject(texDens[deviceid]));
        memmgr.cudaFreeArray(d_dens_array[deviceid]);

        checkCudaErrors(cudaDestroyTextureObject(texSpectrum[deviceid]));
        memmgr.cudaFreeArray(d_spectrum_array[deviceid]);

        checkCudaErrors(cudaDestroyTextureObject(texKern[deviceid]));
        memmgr.cudaFreeArray(d_kern_array[deviceid]);

        if (debugwrite) {
            memmgr.cudaFree(device_dose_accum[deviceid]);
        }
    }

    return true;
}
