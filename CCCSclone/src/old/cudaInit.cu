#include "cudaInit.h"
#include "kernel.h"
#include <cuda_runtime.h>
#include <vector>

cudaArray_t old::d_kern_array;
cudaArray_t old::d_spectrum_array;
cudaArray_t old::d_dens_array;

cudaTextureObject_t old::texSpectrum;
cudaTextureObject_t old::texKern;
cudaTextureObject_t old::texDens;

old::DEVICE_CONV_DATA old::device_data;
float* old::device_dose_accm;

void old::makeTexObject(cudaTextureObject_t *texobj, cudaArray *res, int dims,
    cudaTextureAddressMode addressmode, cudaTextureFilterMode filtermode,
    bool normalizedCoords, 
    cudaTextureReadMode readmode,
    cudaResourceType restype)
{
    cudaChannelFormatDesc channelDesc = 
        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = res;

    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    for (int i=0; i<dims; i++)
        texDesc.addressMode[i] = addressmode;
    texDesc.filterMode = filtermode;
    texDesc.readMode = readmode;
    texDesc.normalizedCoords = normalizedCoords;

    cudaCreateTextureObject(texobj, &resDesc, &texDesc, 0);
}

void old::makeSurfObject(cudaSurfaceObject_t *surfobj, cudaArray *res)
{
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    cudaCreateSurfaceObject(surfobj, &resDesc);
}

int old::initCudaConstandTex(
    SHM_DATA     *datavols,        // host-side data arrays
    MONO_KERNELS *mono,            // spectrum data
    CONSTANTS    *constants,       // calculation paramaters/information
    int          nrays,            // total number of rays
    bool         verbose,    // output switch
    bool         timing,     // output switch
    bool         debug_write // write debug data to file
)
{
    cudaChannelFormatDesc floatChannelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    checkCudaErrors(cudaMallocArray(&d_kern_array, &floatChannelDesc, constants->nradii, constants->ntheta));
    checkCudaErrors(cudaMemcpy2DToArray(d_kern_array, 0, 0, datavols->kernel_array.data(), 
        constants->nradii*sizeof(float), constants->nradii*sizeof(float), constants->ntheta, cudaMemcpyHostToDevice));
    makeTexObject(&texKern, d_kern_array, 2, cudaAddressModeClamp, cudaFilterModeLinear);

    // spectrum data needed for TERMA calculation & normalization
    std::vector<float> spectrum(4 * mono->nkernels);
    // float kermac0 = 0.0f, terma0 = 0.0f;
    for (int e=0; e<mono->nkernels; e++)
    {
        // // calculate T and Kc at zero depth for use in beam hardening correction for terma
        // // see notes in Siddon implementation
        // // zero-depth terma/kerma used in correction of terma for polyenergetic kernel at depth
        // kermac0 += mono->fluence[e]*mono->energy[e]*mono->mu_en[e];
        // terma0 += mono->fluence[e]*mono->energy[e]*mono->mu[e];

        spectrum[e                   ] = mono->fluence[e];
        spectrum[e +   mono->nkernels] = mono->energy[e];
        spectrum[e + 2*mono->nkernels] = mono->mu_en[e];
        spectrum[e + 3*mono->nkernels] = mono->mu[e];
    }
    constants->beamhard_correct = 1.0f; // XXX

    // bind spectrum data to texture memory (with nearest-style fetching - no interpolation)
    // allocated as (nkernels x 4) dim matrix where:
    //   -- (:, 1): fluence
    //   -- (:, 2): energy
    //   -- (:, 3): mu_en
    //   -- (:, 4): mu
    // used in cudaSiddon for terma volume calculation
    checkCudaErrors(cudaMallocArray(&d_spectrum_array, &floatChannelDesc, mono->nkernels, 4));
    checkCudaErrors(cudaMemcpy2DToArray(d_spectrum_array, 0, 0, spectrum.data(), 
        mono->nkernels*sizeof(float), mono->nkernels*sizeof(float), 4, cudaMemcpyHostToDevice));
    makeTexObject(&texSpectrum, d_spectrum_array, 2, cudaAddressModeClamp, cudaFilterModePoint);
    
    // copy mean radii of kernel data to constant memory (first convert radial bounds to mean radii)
    std::vector<float> mean_radii(constants->nradii);
    mean_radii[0] = 0.5 * datavols->radial_boundary[0];
    for (int rr=1; rr<constants->nradii; rr++)
        mean_radii[rr] = 0.5*(datavols->radial_boundary[rr]+datavols->radial_boundary[rr-1]);
    kern_radii_init(mean_radii);

    // Bind denisty to 3D Texture Array
    cudaExtent volExtent = make_cudaExtent(
        constants->size.x,
        constants->size.y,
        constants->size.z
        );
    checkCudaErrors(cudaMalloc3DArray(&d_dens_array, &floatChannelDesc, volExtent));
    cudaMemcpy3DParms CopyParams = {0};
    CopyParams.srcPtr   = make_cudaPitchedPtr(datavols->density.data(), 
        volExtent.width*sizeof(float), volExtent.width, volExtent.height);
    CopyParams.dstArray = d_dens_array;
    CopyParams.extent   = volExtent;
    CopyParams.kind     = cudaMemcpyHostToDevice;
    cudaMemcpy3DAsync(&CopyParams);
    makeTexObject(&texDens, d_dens_array, 3, cudaAddressModeBorder, cudaFilterModeLinear);

    // Find extent of REV data volume
    // Currently this is set to an isotropic box with side length == maximum diagnal length of calc-bbox
    // This will depend on the rotation angles for beam and convolution rays
    std::cout << "Full Convolution Memory Dimensions: " << constants->max_rev_size.x 
        << " x " << constants->max_rev_size.y << " x " << constants->max_rev_size.z << std::endl;
    // allocate intermidate data volume
    int revSize = constants->max_rev_size.x * constants->max_rev_size.y * 
        constants->max_rev_size.z * sizeof(float);
    checkCudaErrors(cudaMalloc((void**)(&(device_data.revDens)), revSize));
    checkCudaErrors(cudaMalloc((void**)(&(device_data.revTerma)), revSize));

    // // Cuda Array for terma texture object fetching
    // cudaExtent calc_bbox_extent = make_cudaExtent(constants->calc_bbox_size.x,
    //     constants->calc_bbox_size.y, constants->calc_bbox_size.z);
    // checkCudaErrors(cudaMalloc3DArray(&(device_data.term_Array), &floatChannelDesc, calc_bbox_extent));
    // makeTexObject(&(device_data.texTerma), 
    //     device_data.term_Array, 3, cudaAddressModeBorder, cudaFilterModeLinear);

    // Allocate 3D CUDA array to bind with surface and texture objects
    floatChannelDesc = cudaCreateChannelDesc<float>();
    cudaExtent dataSize;
    dataSize.width = constants->max_rev_size.x;
    dataSize.height = constants->max_rev_size.y;
    dataSize.depth = constants->max_rev_size.z;
    checkCudaErrors(cudaMalloc3DArray(&(device_data.dose_Array), 
        &floatChannelDesc, dataSize, cudaArraySurfaceLoadStore));
    cudaResourceDesc surfRes;
    memset(&surfRes, 0, sizeof(cudaResourceDesc));
    surfRes.resType = cudaResourceTypeArray;
    surfRes.res.array.array = device_data.dose_Array;
    checkCudaErrors(cudaCreateSurfaceObject(&(device_data.surfDose), &surfRes));
    cudaResourceDesc texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array = device_data.dose_Array;
    cudaTextureDesc texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));
    texDescr.normalizedCoords = false;
    texDescr.filterMode = cudaFilterModeLinear;
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.addressMode[1] = cudaAddressModeBorder;
    texDescr.addressMode[2] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;
    checkCudaErrors(cudaCreateTextureObject(&(device_data.texDose), &texRes, &texDescr, nullptr));
    return 1;
}

int old::freeCudaTexture()
{
    checkCudaErrors(cudaFree(device_data.revDens));
    checkCudaErrors(cudaFree(device_data.revTerma));
    checkCudaErrors(cudaDestroySurfaceObject(device_data.surfDose));
    checkCudaErrors(cudaDestroyTextureObject(device_data.texDose));
    checkCudaErrors(cudaDestroyTextureObject(device_data.texTerma));
    checkCudaErrors(cudaFreeArray(device_data.dose_Array));
    checkCudaErrors(cudaFreeArray(device_data.term_Array));

    checkCudaErrors(cudaDestroyTextureObject(texDens));
    checkCudaErrors(cudaFreeArray(d_dens_array));
    checkCudaErrors(cudaDestroyTextureObject(texSpectrum));
    checkCudaErrors(cudaFreeArray(d_spectrum_array));
    checkCudaErrors(cudaDestroyTextureObject(texKern));
    checkCudaErrors(cudaFreeArray(d_kern_array));
    return 0;
}