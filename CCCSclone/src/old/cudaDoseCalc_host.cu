#include "cudaDoseCalc.h"
#include "cudaInit.h"
#include "kernel.h"
#include "brain_defs.h"
#include "configure.h"
#include "binary_io.h"
#include "geometry.h"
#include "debugLog.h"

#include <iostream>
#include <string>
#include <iomanip>
#include "cuda_runtime.h"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

int old::radconvolveTexture (
    MONO_KERNELS        *mono,
    CONSTANTS           *constants,
    std::vector<BEAM>&  beams,
    int                 nrays,
    RES_LOG& result)
{
    float3 rev_voxelsize = {constants->rev_latspacing, 
        constants->rev_longspacing, constants->rev_latspacing};

    // Set up arrays for dynamic GPU resource allocation
    // depending on dimensions of BEV data per convolution ray
    dim3 tileBlock(TILE_DIM_X, TILE_DIM_Y, 1);
    std::vector<dim3> rayGrid(nrays);
    std::vector<dim3> conBlock(nrays);
    std::vector<dim3> conGrid(nrays);
    std::vector<uint> memsize(nrays);

    // Calculate cuda execution block/grid sizes
    uint dataSize = constants->nvoxels();
    uint calcDataSize = constants->bbox_nvoxels();

    uint3 max_actual_rev_size = {0};
    for (int dc=0; dc<beams.size(); dc++)
    {
        BEAM& this_beam = beams[dc];
        std::cout << "Starting: | beam#" << dc << std::endl;

        // 1. initialize fluence map. Here, we enable all beamlets
        float* d_fluence_map;
        int fluenceSize = this_beam.fmap_size.x * this_beam.fmap_size.y;
        std::vector<float> h_fluence_map(fluenceSize, 1.);
        checkCudaErrors(cudaMalloc((void**)(&d_fluence_map), fluenceSize * sizeof(float)));
        checkCudaErrors(cudaMemcpy(d_fluence_map, h_fluence_map.data(), 
            fluenceSize*sizeof(float), cudaMemcpyHostToDevice));

        PILLAR_GRID hPG{};
        hPG.numBeamlets = this_beam.fmap_size.x * this_beam.fmap_size.y;
        hPG.beamletIdx = std::vector<int>(hPG.numBeamlets, 0);
        hPG.pillarStartCoords = std::vector<float3>(hPG.numBeamlets);
        hPG.beamletAngles = std::vector<float2>(hPG.numBeamlets);
        hPG.beamletIsocenters = std::vector<float3>(hPG.numBeamlets);

        for (int i=0; i<hPG.numBeamlets; i++)
            hPG.beamletIdx[i] = i;
        
        // 2. get the central axis limit (anchors) of beamlets
        std::vector<float3> beamlet_start(hPG.numBeamlets);
        std::vector<float3> beamlet_end(hPG.numBeamlets);
        std::vector<float> beamlet_length(hPG.numBeamlets);
        float max_beamlet_length = 1.0f;
        hPG.max_beamlet_size = float2{0.f, 0.f};
        for (int i=0; i<hPG.numBeamlets; i++)
        {
            calcBeamletAnchors(
                beamlet_start[i], beamlet_end[i], hPG.beamletAngles[i], hPG.beamletIsocenters[i],
                this_beam.source, this_beam.isocenter,
                static_cast<uint>(hPG.beamletIdx[i]),
                this_beam.beamlet_size, this_beam.fmap_size,
                constants->voxel,
                constants->start,
                constants->calc_bbox_start,
                constants->calc_bbox_size,
                this_beam.azimuth, this_beam.zenith, this_beam.coll
            );
            beamlet_length[i] = length(beamlet_end[i] - beamlet_start[i]);
            max_beamlet_length = std::max(max_beamlet_length, beamlet_length[i]);
            float2 beamlet_diverge_size = this_beam.beamlet_size * 
                length(beamlet_end[i]-this_beam.source) /
                length(hPG.beamletIsocenters[i] - this_beam.source);
            hPG.max_beamlet_size.x = max(hPG.max_beamlet_size.x, beamlet_diverge_size.x);
            hPG.max_beamlet_size.y = max(hPG.max_beamlet_size.y, beamlet_diverge_size.y);
        }

        // compute pillar size and ensure that it is an integer multiple of rev_voxelsize
        float psize_long = max_beamlet_length + 2. * constants->kernel_extent + hPG.wallThickness*rev_voxelsize.y;
        float2 psize_trans = hPG.max_beamlet_size + 2. * constants->kernel_extent + hPG.wallThickness*rev_voxelsize.x;
        float3 expand = rev_voxelsize - make_float3(
            fmodf(psize_trans.x, rev_voxelsize.x),
            fmodf(psize_long,    rev_voxelsize.y),
            fmodf(psize_trans.y, rev_voxelsize.z));
        hPG.pillarSize = make_float3(psize_trans.x, psize_long, psize_trans.y) + 
            expand + 2.*float(hPG.pillarBuffer)*rev_voxelsize;
        hPG.pillarDims = make_int3(hPG.pillarSize / rev_voxelsize);

        // compute pillar limits to use in geometric transformations
        for (int i=0; i<hPG.numBeamlets; i++)
        {
            float3 g_offset = make_float3(
                -0.5f * (hPG.pillarSize.x + hPG.wallThickness*rev_voxelsize.x),
                -0.5f * (hPG.pillarSize.y - beamlet_length[i] + hPG.wallThickness*rev_voxelsize.y),
                -0.5f * (hPG.pillarSize.z + hPG.wallThickness*rev_voxelsize.z) );
            // TODO fix rotation pivot to src here (may be unnecessary)
            g_offset = inverseRotateBeamAtOriginRHS(
                    inverseRotateBeamletAtOriginRHS(g_offset, hPG.beamletAngles[i].x, hPG.beamletAngles[i].y),
                        this_beam.azimuth, this_beam.zenith, this_beam.coll);
            hPG.pillarStartCoords[i] = beamlet_start[i] + g_offset;
        }

        // 3. Build the packed BEV allocation map
        hPG.numPillars.x = static_cast<int>(std::ceil(std::sqrt(float(hPG.numBeamlets))));
        hPG.numPillars.y = static_cast<int>(std::ceil(float(hPG.numBeamlets) / hPG.numPillars.x));

        // add another wall at the end of the array in each lateral dimension
        hPG.gridDims = make_int3(
            hPG.numPillars.x * hPG.pillarDims.x + hPG.wallThickness,
            hPG.pillarDims.y + hPG.wallThickness,
            hPG.numPillars.y * hPG.pillarDims.z + hPG.wallThickness
        );
        hPG.gridSize = make_float3(hPG.gridDims) * rev_voxelsize;

        int*    dpg_beamletIdx;
        float3* dpg_pillarStartCoords;
        float2* dpg_beamletAngles;
        float3* dpg_beamletIsocenters;
        checkCudaErrors(cudaMalloc((void**)(&dpg_beamletIdx), hPG.numBeamlets*sizeof(int)));
        checkCudaErrors(cudaMalloc((void**)(&dpg_pillarStartCoords), hPG.numBeamlets*sizeof(float3)));
        checkCudaErrors(cudaMalloc((void**)(&dpg_beamletAngles), hPG.numBeamlets*sizeof(float2)));
        checkCudaErrors(cudaMalloc((void**)(&dpg_beamletIsocenters), hPG.numBeamlets*sizeof(float3)));
        checkCudaErrors(cudaMalloc((void**)(&(device_data.dose)), hPG.pillar_grid_nvoxels()*sizeof(float)));

        // Move pillar_grid data to device
        checkCudaErrors(cudaMemcpy(dpg_beamletIdx, hPG.beamletIdx.data(), hPG.numBeamlets*sizeof(int), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(dpg_pillarStartCoords, hPG.pillarStartCoords.data(), hPG.numBeamlets*sizeof(float3), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(dpg_beamletAngles, hPG.beamletAngles.data(), hPG.numBeamlets*sizeof(float2), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(dpg_beamletIsocenters, hPG.beamletIsocenters.data(), hPG.numBeamlets*sizeof(float3), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemset(device_data.dose, 0., hPG.pillar_grid_nvoxels()*sizeof(float)));

        dim3 packedGrid = dim3(
            static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.y) / tileBlock.x)),
            static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.z) / tileBlock.y)),
            static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.x) / tileBlock.z))
        );

        dim3 rayBlock(1, 32, 3);
        std::vector<REV_DATA> rev(nrays);

        // compute
        radconvolvePrep(constants, &hPG,
            nrays, rev, this_beam, rayGrid, rayBlock,
            conGrid, conBlock, memsize, max_actual_rev_size);

        radconvolveCompute(
            mono, constants, nrays, this_beam,
            hPG, dpg_beamletIdx, dpg_beamletAngles,
            dpg_pillarStartCoords, dpg_beamletIsocenters,
            rayGrid, rayBlock,
            rev, d_fluence_map,
            conGrid, conBlock,
            memsize, packedGrid, tileBlock,
            dc, result[dc]
        );

        // clean up
        checkCudaErrors(cudaFree(device_data.dose));
        checkCudaErrors(cudaFree(dpg_beamletIsocenters));
        checkCudaErrors(cudaFree(dpg_beamletAngles));
        checkCudaErrors(cudaFree(dpg_pillarStartCoords));
        checkCudaErrors(cudaFree(dpg_beamletIdx));
        checkCudaErrors(cudaFree(d_fluence_map));
    }
    return 0;
}

int old::radconvolvePrep(CONSTANTS* constants, PILLAR_GRID* hPG, 
    int nrays, std::vector<REV_DATA>& rev, BEAM& this_beam, std::vector<dim3>& rayGrid, dim3& rayBlock,
    std::vector<dim3>& conGrid, std::vector<dim3>& conBlock, std::vector<uint>& memsize,
    uint3& max_actual_rev_size)
{
    for (int rr=0; rr<nrays; rr++)
    {
        float theta = constants->get_theta_from_index(rr);
        float phi = constants->get_phi_from_index(rr);

        findREV(&rev[rr],
            constants,
            &this_beam,
            make_float3(0.f),
            hPG->gridSize,
            theta,
            phi);
        
        // PAY ATTENTION: rev[rr].size in the order of YZX
        // rayGrid[rr] in the order of XYZ
        rayGrid[rr] = make_uint3(rev[rr].size.z,
            static_cast<uint>(std::ceil(static_cast<float>(rev[rr].size.x) / rayBlock.y)),
            static_cast<uint>(std::ceil(static_cast<float>(rev[rr].size.y) / rayBlock.z)));
        
        // rev[rr].size in the order of YZX
        if (rev[rr].size.x > 1024)
            throw std::runtime_error("size of rev (x-dim) exceeds maximum #threads allowed. Contact dev to resolve.");
        
        // buffer beyond used rev block, where convolve kernel 
        // will also write zero before proceeding
        // conBlock in the order of YZX
        int overwrite = 2;
        conBlock[rr].x = rev[rr].size.x + overwrite;
        conGrid[rr].y  = rev[rr].size.y + overwrite;
        conGrid[rr].z  = rev[rr].size.z + overwrite;
        conBlock[rr].y = conBlock[rr].z = conGrid[rr].x = 1;
        memsize[rr] = 2 * (conBlock[rr].x - overwrite) * sizeof(float);

        max_actual_rev_size.x = std::max(max_actual_rev_size.x, rev[rr].size.x);
        max_actual_rev_size.y = std::max(max_actual_rev_size.y, rev[rr].size.y);
        max_actual_rev_size.z = std::max(max_actual_rev_size.z ,rev[rr].size.z);
    }
    return 0;
}

int old::radconvolveCompute(
    MONO_KERNELS *mono, CONSTANTS* constants, int nrays, BEAM& this_beam,
    PILLAR_GRID& hPG, int* dpg_beamletIdx, float2* dpg_beamletAngles,
    float3* dpg_pillarStartCoords, float3* dpg_beamletIsocenters,
    const std::vector<dim3>& rayGrid, const dim3& rayBlock,
    const std::vector<REV_DATA>& rev, float* d_fluence_map,
    const std::vector<dim3>& conGrid, const std::vector<dim3>& conBlock, 
    const std::vector<uint>& memsize, const dim3& packedGrid, const dim3& tileBlock,
    int dc, BEAM_LOG& beamResult)
{
    float3 rev_voxelsize = {constants->rev_latspacing, constants->rev_longspacing, constants->rev_latspacing};
    fs::path debugDir(Paths::Instance()->debug_dir());
    for (int ray_idx=0; ray_idx<nrays; ray_idx++)
    {
        // get kernel rotation angles once instead of for each kernel thread launch
        int kern_wt_idx = constants->get_kernel_theta_index(ray_idx);
        float kern_theta = constants->get_theta_from_index(ray_idx);
        float kern_phi = constants->get_phi_from_index(ray_idx);

        // calculate rev Terma and sample density
        // rayGrid[rr] in the order of XYZ
        int raySharedMem = rayBlock.x * rayBlock.y * rayBlock.z * sizeof(float);

        float3* g_coords_log = nullptr;
        if (constants->debugREVTerma && dc==0 && ray_idx==0)
        {
            #include "debug_raytrace.cpp.in"
            int volume = constants->max_rev_size.x * 
                constants->max_rev_size.y * constants->max_rev_size.z;
            checkCudaErrors(cudaMalloc((void**)(&g_coords_log), volume*sizeof(float3)));
        }

        cudaBeamletRaytrace<<<rayGrid[ray_idx], rayBlock, raySharedMem>>>(
            device_data.revDens,
            device_data.revTerma,
            this_beam.source,
            this_beam.beamlet_size,
            this_beam.azimuth,
            this_beam.zenith,
            this_beam.coll,
            make_float3(hPG.gridDims),
            make_float3(hPG.pillarDims),
            hPG.numBeamlets,
            hPG.wallThickness,
            hPG.numPillars,
            dpg_beamletIdx,
            dpg_beamletAngles,
            dpg_pillarStartCoords,
            dpg_beamletIsocenters,
            kern_theta, kern_phi,
            rev[ray_idx].min_coords,
            rev[ray_idx].size,
            constants->max_rev_size,
            constants->start,
            constants->voxel,
            rev_voxelsize,
            make_float3(constants->calc_bbox_start),
            texDens,
            g_coords_log,

            // raytrcing/terma args
            d_fluence_map,
            make_float3(constants->size),
            constants->calc_bbox_size,
            constants->beamhard_correct,
            this_beam.direction,
            this_beam.isocenter,
            this_beam.sad,
            this_beam.fmap_size,
            mono->nkernels,
            texSpectrum
        );

        float* d_debugProbe=nullptr;
        if (constants->debugREVTerma && dc==0 && ray_idx==0)
        {
            if (write_g_coords_log(g_coords_log, constants, debugDir))
                return 1;
            checkCudaErrors(cudaFree(g_coords_log));
            if (writeREVTerma(device_data.revTerma, 
                device_data.revDens, constants, debugDir))
                return 1;
            #include "debug_convolve.cpp.in"
        }
        if (constants->debugREVDose && dc==0 && ray_idx==0)
        {
            int volume = constants->max_rev_size.x * 
                constants->max_rev_size.y * constants->max_rev_size.z;
            checkCudaErrors(cudaMalloc((void**)(&d_debugProbe), volume*sizeof(float)));
            checkCudaErrors(cudaMemset(d_debugProbe, 0, volume*sizeof(float)));
        }

        // perform dose calculation (CCCS) w/ heterogeneity correction in REV volume
        PackRowConvolve<<<conGrid[ray_idx], conBlock[ray_idx], memsize[ray_idx]>>> (
            device_data.revDens,
            device_data.revTerma,
            device_data.surfDose,
            (float)kern_wt_idx,
            rev[ray_idx].size,
            constants->max_rev_size,
            constants->rev_longspacing,
            constants->nradii,
            constants->ntheta,
            constants->nphi,
            texKern,
            d_debugProbe
        );

        if (constants->debugREVDose && dc==0 && ray_idx==0)
        {
            if (writeREVDose(device_data.texDose, constants, debugDir))
                return 1;
            if (writeREVDebug(d_debugProbe, constants, debugDir))
                return 1;
            if (writeREVSurf(device_data.surfDose, constants, debugDir))
                return 1;
            checkCudaErrors(cudaFree(d_debugProbe));
        }

        // transform packed REV dose coefficients from the previous convolution back to the BEV system then
        // perform element-by-element sum, accumulating over all convolution directions
        PackedREVtoBEVdose <<<packedGrid, tileBlock, 0>>> (
            device_data.dose,             // beamlet-packed dose array in BEV orientation
            device_data.texDose,          // packed dose array embedded in REV bounding box
            kern_theta, kern_phi,         // convolution direction
            rev[ray_idx].min_coords,      // REV volume limit coords in XYZ coord system
            rev_voxelsize,
            hPG.gridDims
        );

        // for debug purposes
        if (constants->debugREVDose)
            std::cout << "ray index: " << ray_idx << std::endl;
    }

    if (constants->debugBEVDose && dc==0)
        if (writeBEVDose(device_data.dose, hPG, debugDir))
            return 1;
    
    // if (! constants->logPatientDose)
    //     return 0;

    // unpack pillars from BEV storage
    // copy output of PackedREVtoBEVdose to cudaArray and attach texture object
    cudaChannelFormatDesc floatChannelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaExtent packedArrayExtent = make_cudaExtent(hPG.gridDims.y, hPG.gridDims.z, hPG.gridDims.x);
    cudaArray_t PackedBEVdose_Array;
    checkCudaErrors(cudaMalloc3DArray(&PackedBEVdose_Array, &floatChannelDesc, packedArrayExtent));
    // copy to cudaArray
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr((void*)device_data.dose, 
        packedArrayExtent.width*sizeof(float), packedArrayExtent.width, 
        packedArrayExtent.height);
    copyParams.dstArray = PackedBEVdose_Array;
    copyParams.extent = packedArrayExtent;
    copyParams.kind = cudaMemcpyDeviceToDevice;
    checkCudaErrors(cudaMemcpy3DAsync(&copyParams));
    // attach to texture object
    cudaTextureObject_t texPackedBEVDose;
    makeTexObject(&texPackedBEVDose, PackedBEVdose_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);

    // determine kernel launch params
    dim3 unpackBlock = tileBlock;
    dim3 unpackGrid = dim3(
        static_cast<uint>(std::ceil(static_cast<float>(constants->calc_bbox_size.x) / unpackBlock.x)),
        static_cast<uint>(std::ceil(static_cast<float>(constants->calc_bbox_size.y) / unpackBlock.y)),
        static_cast<uint>(std::ceil(static_cast<float>(constants->calc_bbox_size.z) / unpackBlock.z))
    );

    int calcDataSize = constants->bbox_nvoxels();
    std::get<1>(beamResult).resize(hPG.numBeamlets);
    float* d_unpacked_dose;
    checkCudaErrors(cudaMalloc((void**)(&d_unpacked_dose), calcDataSize*sizeof(float)));
    for (int i=0; i<hPG.numBeamlets; i++)
    {
        std::get<1>(beamResult)[i].resize(calcDataSize);
        checkCudaErrors(cudaMemset(d_unpacked_dose, 0., calcDataSize*sizeof(float)));
        UnpackBEVDosePillar <<<unpackGrid, unpackBlock>>> (
            d_unpacked_dose,
            texPackedBEVDose,
            this_beam.sad,
            this_beam.source,
            this_beam.azimuth,
            this_beam.zenith,
            this_beam.coll,
            constants->start,
            constants->voxel,
            make_float3(constants->calc_bbox_start),
            constants->calc_bbox_size,
            rev_voxelsize,
            i,
            i % hPG.numPillars.x,
            i / hPG.numPillars.x,
            hPG.pillarDims,
            hPG.wallThickness,
            hPG.pillarBuffer,
            hPG.pillarStartCoords[i],
            hPG.beamletAngles[i]
        );
        checkCudaErrors(cudaMemcpy(std::get<1>(beamResult)[i].data(), d_unpacked_dose, 
            calcDataSize*sizeof(float), cudaMemcpyDeviceToHost));
    }
    std::get<0>(beamResult) = hPG;
    checkCudaErrors(cudaFree(d_unpacked_dose));
    checkCudaErrors(cudaDestroyTextureObject(texPackedBEVDose));
    checkCudaErrors(cudaFreeArray(PackedBEVdose_Array));

    return 0;
}