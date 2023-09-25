#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <cmath>
#include <cassert>
#include <utility>
#include <vector>
#include <exception.h>

#if WITH_OPENCV2
#include <opencv2/core/core.hpp>
#else
#include <opencv2/imgproc.hpp>
#endif

#include "../include/debug.h"
#include "cudaMemFunctions.cuh"
#include "kernels.cuh"
#include "DoseCalcIO/volume.h"
#include "DoseCalcIO/binary_io.h"
#include "DoseCalcIO/roi.h"
#include "CudaUtilities/vector_ops.h"
#include "CudaUtilities/make_tex_surf.cuh"
#include "./debug.h"
#include "CudaUtilities/array_sum.cuh"

extern bool debugwrite;

// take HU data volume as input
// convert to density, resize and resample according to densize, with isotropic resolution
int read_dicom(FloatVolume& dens, FloatVolume& ctdata, float iso_voxel, CTLUT* ctlut, bool verbose)
{
    cudaCreateTexIso( dens, ctdata, iso_voxel, ctlut );
    if (verbose) {
        printf("RESAMPLED DENSITY ATTRIBUTES\n");
        printf("----------------------------\n");
        printf("Density dimensions:         %d x %d x %d\n",dens.size.x,dens.size.y,dens.size.z);
        printf("Density voxel size [cm]:    %3.2f x %3.2f x %3.2f\n",dens.voxsize.x,dens.voxsize.y,dens.voxsize.z);
        printf("Density volume origin [cm]: %5.2f , %5.2f , %5.2f\n\n",dens.start.x,dens.start.y,dens.start.z);
    }
    return 1;
}

// convert digital phantoms or cti phantoms from HU to density
void cudaHU2dens( FloatVolume& dens, FloatVolume& origCT, int gS, bool verbose )
{
	// data sizes for allocating GPU memory
    const int shortSize = origCT.mem_size();
    const int floatSize = dens.mem_size();

	// allocate GPU array to hold MU volume
    float *d_origCT;
    checkCudaErrors(cudaMalloc( (void**) &d_origCT, shortSize ));
    checkCudaErrors(cudaMemcpy( d_origCT, origCT.data(), shortSize, cudaMemcpyHostToDevice ) );

	// allocate GPU array to hold density volume
    float *dens_ray;
    checkCudaErrors(cudaMalloc( (void**) &dens_ray, floatSize ));

    // 512 threads per block
    dim3 block(512);
    // determine grid size (number of blocks) according to total data size / threads per block
    unsigned int sizer = dens.nvoxels();
    unsigned int gzero = sizer / block.x;
    if (sizer % block.x > 0) gzero++; // check for data size != multiple of 512
    dim3 grid(gzero);
    if (verbose)
        printf("\n gridSize: %d", grid.x);

    // cuda kernel implemented in "cudaCalcDeff.cuh"
    deviceHU2dens<<<grid,block>>>( d_origCT, dens_ray, dens.size );
    getLastCudaError("Kernel execution failed");

	// copy results back to CPU
    checkCudaErrors( cudaMemcpy( dens.data(), dens_ray, floatSize, cudaMemcpyDeviceToHost ) );

	// free GPU arrays
    cudaFree( d_origCT );
    cudaFree( dens_ray );
}

// resize and resample data into an isotropic volume with dimensions of densize x densize x densize voxels
// isotropic volume resolution automatically determined by densize dimensions
void cudaCreateTexIso( FloatVolume& dens, FloatVolume& data, float iso_voxel, CTLUT* ctlut, bool hu2dens)
{
    /* if ctlut is valid ptr, HU values are converted to density values in kernel: cudaMakeIsotropic */
    // get GPU properties
    cudaDeviceProp devProp;
    cudaGetDeviceProperties( &devProp, 0 );

    /////////////////// Bind Inputs to 3D Texture Arrays //////////////////////////////////////////////
    // bind the HU data volume (data.data) to a 3D texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent imgSize	= make_cudaExtent((float)data.size.x, (float)data.size.y, (float)data.size.z);
    cudaMalloc3DArray(&imgArray, &channelDesc, imgSize);

    cudaMemcpy3DParms imgParams = {0};
    imgParams.srcPtr	    =	make_cudaPitchedPtr((void*)data.data(), imgSize.width*sizeof(float), imgSize.width, imgSize.height);
    imgParams.dstArray	    =	imgArray;
    imgParams.extent	    =	imgSize;
    imgParams.kind		    =	cudaMemcpyHostToDevice;
    cudaMemcpy3D(&imgParams);

    cudaResourceDesc texImgRes;
    memset(&texImgRes, 0, sizeof(cudaResourceDesc));
    texImgRes.resType = cudaResourceTypeArray;
    texImgRes.res.array.array = imgArray;

    // texImg.normalized      	=	false;
    // texImg.filterMode	    =	cudaFilterModeLinear;
    // texImg.addressMode[0]	=	cudaAddressModeBorder;
    // texImg.addressMode[1]	=	cudaAddressModeBorder;
    // texImg.addressMode[2]	=	cudaAddressModeBorder;

    cudaTextureDesc texImgDescr;
    memset(&texImgDescr, 0, sizeof(cudaTextureDesc));
    texImgDescr.normalizedCoords = false;
    texImgDescr.filterMode = cudaFilterModeLinear;
    texImgDescr.addressMode[0] = cudaAddressModeBorder;
    texImgDescr.addressMode[1] = cudaAddressModeBorder;
    texImgDescr.addressMode[2] = cudaAddressModeBorder;

    // cudaBindTextureToArray(texImg, imgArray, channelDesc);
    checkCudaErrors(cudaCreateTextureObject(&texImg, &texImgRes, &texImgDescr, NULL));
    //////////////////////////////////////////////////////////////////////////////////////////////

    // calculate new dimensions
    uint3 iso_size = float2uint_ceil(data.voxsize*(make_float3(data.size)/make_float3(iso_voxel)));
    uint iso_count = product(iso_size);
    uint iso_memsize = iso_count * sizeof(float);

	// allocate GPU array to hold isotropic results
    checkCudaErrors( cudaMalloc( (void**) &iso_matrix, iso_memsize ) );
    checkCudaErrors( cudaMemset( iso_matrix, 0, iso_memsize ) );

    // set number of threads and blocks
    dim3 block( devProp.maxThreadsPerBlock );
    unsigned int sizer = iso_count;

    int3 gridSize;
    gridSize.x = sizer / block.x;
    gridSize.y = 1;
    gridSize.z = 1;
    if ( sizer % block.x > 0 ) gridSize.x++;
    if ( gridSize.x > devProp.maxGridSize[1] )
    {
        gridSize.y = gridSize.x / devProp.maxGridSize[1];
        if ( gridSize.x % devProp.maxGridSize[1] > 0 ) gridSize.y++;
        gridSize.x = devProp.maxGridSize[1];

        if ( gridSize.y > devProp.maxGridSize[1] )
        {
            gridSize.z = gridSize.y / devProp.maxGridSize[1];
            if ( gridSize.y % devProp.maxGridSize[1] > 0 ) gridSize.z++;
            gridSize.y = devProp.maxGridSize[1];
        }
    }
    dim3 grid( gridSize.x, gridSize.y, gridSize.z );

    if (ctlut != nullptr) {
        uint lutsize = ctlut->points.size();

        // make linterp texture for use in kernel
        float* h_hunits = new float[lutsize];
        float* h_massdens = new float[lutsize];
        for (int ii=0; ii<lutsize; ii++) {
            h_hunits[ii] = ctlut->points[ii].hunits;
            h_massdens[ii] = ctlut->points[ii].massdens;
        }
        float* d_hunits;
        float* d_massdens;
        checkCudaErrors( cudaMalloc((void**)&d_hunits, lutsize*sizeof(*d_hunits)) );
        checkCudaErrors( cudaMemcpy(d_hunits, h_hunits, lutsize*sizeof(*h_hunits), cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMalloc((void**)&d_massdens, lutsize*sizeof(*d_massdens)) );
        checkCudaErrors( cudaMemcpy(d_massdens, h_massdens, lutsize*sizeof(*h_massdens), cudaMemcpyHostToDevice) );
        int isoShared = 2*sizeof(float)*lutsize;

        // run cuda kernel to create isotropic density volume with axes voxel count: densize
        cudaMakeIsotropicWithLUT<<< grid, block, isoShared>>>( iso_matrix, data.voxsize, iso_voxel, iso_size, 
            d_hunits, d_massdens, lutsize, texImg);
        getLastCudaError("cudaMakeIsotropic()");

        delete [] h_hunits;
        delete [] h_massdens;
        checkCudaErrors(cudaFree(d_hunits));
        checkCudaErrors(cudaFree(d_massdens));
    } else {
        // run cuda kernel to create isotropic density volume with axes voxel count: densize
        cudaMakeIsotropic<<< grid, block>>>( iso_matrix, data.voxsize, iso_voxel, iso_size, texImg);
        if (hu2dens) {
            deviceHU2dens<<<grid, block>>>(iso_matrix, iso_matrix, iso_size);
        }
        getLastCudaError("cudaMakeIsotropic()");
    }

    // checkCudaErrors(cudaUnbindTexture(texImg));
    // checkCudaErrors(cudaFreeArray(imgArray));
    if (texImg)
        checkCudaErrors(cudaDestroyTextureObject(texImg));
    
    if (imgArray)
        checkCudaErrors(cudaFreeArray(imgArray));

    // copy array properties
    dens.size = iso_size;
    dens.voxsize = make_float3(iso_voxel);
    dens.start = data.start;

    //allocate density array and copy results from GPU to CPU
    dens.set_data(dens.nvoxels());
    checkCudaErrors(cudaMemcpy(dens.data(),iso_matrix,iso_memsize,cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(iso_matrix));
}

Volume<uint8_t> generateContourMask(StructureSet& contour, FrameOfReference ctframe, FrameOfReference densframe, void* texRay ) {
  // iterate over axial slices, converting contour point coords list to closed polygon
  // Construct volume in original dicom coordsys from rasterized slices
  // resample dicom coordsys to dose coordsys
  Volume<uint8_t> mask = Volume<uint8_t>(ctframe);

  // for each slice in contour point set, determine where in total volume it exists
  unsigned long mark = 0;
  for (int c=0; c<contour.sub_cntr_count; c++) {
    float cntr_z = 0.1f*contour.points[3*mark+2];
    int z_offset = round((cntr_z - mask.start.z) / mask.voxsize.z);

    // use opencv to draw polygon
    // make sure to += to existing slice data to handle multiple contours in same slice
    // binary thresholding will be done later
    cv::Mat mat = cv::Mat::zeros(mask.size.y, mask.size.x, CV_8UC1);
    cv::Point* pts = new cv::Point[contour.sub_cntr_points_count[c]];

    // generate point set and fill polygon
    for (int p=0; p<contour.sub_cntr_points_count[c]; p++) {
      int idx_x = floor((0.1f*contour.points[3*(mark+p)    ] - mask.start.x) / mask.voxsize.x);
      int idx_y = floor((0.1f*contour.points[3*(mark+p) + 1] - mask.start.y) / mask.voxsize.y);
      pts[p] = cv::Point{idx_x, idx_y};
    }

    // render polygon
    cv::Point* pts_array[1] = { pts };
    const int npts[1] = { (int)contour.sub_cntr_points_count[c] };
    cv::fillPoly(mat, (const cv::Point**)pts_array, npts, 1, cv::Scalar(1));
    delete[] pts;

    // copy polygon to mask volume
    long int memoffset = sizeof(uint8_t)*mask.size.x*mask.size.y*z_offset;
    for (int r=0; r<mat.rows; r++) {
      for (int c=0; c<mat.cols; c++) {
        int ptr = memoffset+r*mat.cols+c;
        mask._vect[ptr] = (uint8_t)((mat.at<uint8_t>(r, c)+mask._vect[ptr])>0);
      }
    }
    mark += contour.sub_cntr_points_count[c];
  }


  // bind the mask to a texture for resampling
  /////////////////// Bind Inputs to 3D Texture Arrays //////////////////////////////////////////////
  cudaArray* arrayMask;
  cudaTextureObject_t texMask;
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<uint8_t>();
  cudaExtent extent  = make_cudaExtent(mask.size.x, mask.size.y, mask.size.z);
  checkCudaErrors( cudaMalloc3DArray(&arrayMask, &channelDesc, extent) );

  cudaMemcpy3DParms parms = {0};
  parms.srcPtr   = make_cudaPitchedPtr((void*)mask.data(), extent.width*sizeof(uint8_t), extent.width, extent.height);
  parms.dstArray = arrayMask;
  parms.extent   = extent;
  parms.kind     = cudaMemcpyHostToDevice;

  cudaMemcpy3D(&parms);
  makeTexObject<cudaArray>(&texMask, arrayMask, 3, cudaAddressModeBorder, cudaFilterModePoint);
  //////////////////////////////////////////////////////////////////////////////////////////////

  uint8_t* d_remask;
  checkCudaErrors( cudaMalloc(&d_remask, sizeof(uint8_t)*densframe.nvoxels()) );
  dim3 rsBlock = dim3{16,16,1};
  dim3 rsGrid = dim3{
    (uint)ceilf(densframe.size.x/rsBlock.x),
    (uint)ceilf(densframe.size.y/rsBlock.y),
    (uint)ceilf(densframe.size.z/rsBlock.z),
  };
  cudaResample<<<rsGrid, rsBlock>>>(
      d_remask, densframe.start, densframe.size, densframe.spacing,
      texMask, ctframe.start,    ctframe.size,   ctframe.spacing
      );
  cudaThreshold<<<rsGrid, rsBlock>>>(
      d_remask, d_remask, densframe.size, (uint8_t)1
      );

  Volume<uint8_t> remask(densframe);
  checkCudaErrors( cudaMemcpy((void*)remask.data(), d_remask, sizeof(uint8_t)*densframe.nvoxels(), cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaFree(d_remask) );
  checkCudaErrors( cudaDestroyTextureObject(texMask) );
  checkCudaErrors( cudaFreeArray(arrayMask) );

  return remask;
}


// projects from beam source through binary contour volume created in "create_contour_volume"
// and outputs a 2D binary mask with resolution equal to fluence map beamlet size (anisotropic)
// perform a raytracing through the volume from the source point of each beam
// and outputs a 2D binary fluence mask defining a conformal field at isocenter
void findFluenceProjection(
        std::vector<float>& fluence_map,
        const Volume<uint8_t>& ptv_mask,
        float3 isocenter,
        float3 source,
        uint2   fmap_size,
        float2 beamlet_size,
        float azimuth,
        float zenith,
        float coll,
        int verbose
        ) {

    // convert mask to float
    std::vector<float> mask_data = std::vector<float>(ptv_mask._vect.begin(), ptv_mask._vect.end());

    // bind the mask to a texture to be used in the creation of the 2D binary fluence masks
    /////////////////// Bind Inputs to 3D Texture Arrays //////////////////////////////////////////////
    cudaArray *rayArray;
    cudaTextureObject_t texRay;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent raySize  = make_cudaExtent(ptv_mask.size.x, ptv_mask.size.y, ptv_mask.size.z);
    cudaMalloc3DArray(&rayArray, &channelDesc, raySize);

    cudaMemcpy3DParms rayParams = {0};
    rayParams.srcPtr        =   make_cudaPitchedPtr((void*)mask_data.data(), raySize.width*sizeof(float), raySize.width, raySize.height);
    rayParams.dstArray      =   rayArray;
    rayParams.extent        =   raySize;
    rayParams.kind          =   cudaMemcpyHostToDevice;

    cudaMemcpy3D(&rayParams);

    makeTexObject<cudaArray>(&texRay, rayArray, 3, cudaAddressModeBorder, cudaFilterModeLinear);
    //////////////////////////////////////////////////////////////////////////////////////////////

    // projection through iso_cntr_matrix to create conformal field map from source
    float *d_fluence;
    int nbixels = fmap_size.x * fmap_size.y;
    checkCudaErrors( cudaMalloc( &d_fluence, nbixels*sizeof(float) ) );
    checkCudaErrors( cudaMemset( d_fluence, 0, nbixels*sizeof(float) ) );

    // bool: should we supersample the raytracing grid?
    bool supersample = true;

    dim3 projBlock;
    dim3 projGrid = dim3{1, 0, 0};
    if (supersample) {
        // use center/corner/edge sampling pattern
        projBlock = dim3{9, 8, 8};
    } else {
        projBlock = dim3{1, 16, 16};
    }
    // create a thread for each pixel in a 2D square fluence map array
    projGrid.y = ceilf((float)fmap_size.x/projBlock.y);
    projGrid.z = ceilf((float)fmap_size.y/projBlock.z);

    if (verbose>=3) {
        printf(" ### projGrid:  %3d x %3d x %3d  |  projBlock:  %3d x %3d x %3d\n",
                projGrid.x, projGrid.y, projGrid.z,
                projBlock.x, projBlock.y, projBlock.z
              );
    }
    // cuda raytracing function implemented in "cudaCalcDeff.cuh"
    cudaRaycast<<<projGrid, projBlock, projBlock.x*projBlock.y*projBlock.z*sizeof(char)>>> (
            d_fluence,
            ptv_mask.start,
            ptv_mask.size,
            ptv_mask.voxsize,
            source,
            isocenter,
            fmap_size,
            beamlet_size,
            azimuth,
            zenith,
            coll,
            texRay
            );
    getLastCudaError("cudaRaycast");

    // copy the results back to the CPU
    checkCudaErrors( cudaDestroyTextureObject(texRay) );
    checkCudaErrors( cudaFreeArray(rayArray) );
    checkCudaErrors( cudaMemcpy( fluence_map.data(), d_fluence, nbixels*sizeof(float), cudaMemcpyDeviceToHost ) );
    checkCudaErrors(cudaFree(d_fluence));
}
