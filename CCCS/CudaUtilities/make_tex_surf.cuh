#ifndef __MAKE_TEX_SURF_CUH__
#define __MAKE_TEX_SURF_CUH__

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cstring>

/* Quick Notes about Texture Objects and Surface Objects
 * For linear textures (where T!=cudaArray), data is read using tex1Dfetch()
 *   and the construction options are:
 *   * cudaTextureAddressMode:
 *     - cudaTextureAddressMode::cudaAddressModeClamp    = 1
 *     - cudaTextureAddressMode::cudaAddressModeBorder   = 4
 *   * cudaTextureFilterMode:
 *     - cudaTextureFilterMode::cudaFilterModePoint      = 0
 *   * cudaTextureReadMode:
 *     - cudaTextureReadMode::cudaReadModeElementType    = 0
 *     - cudaTextureReadMode::cudaReadModeNormalizeFloat = 0
 *   * cannot use normalized coordinate addressing
 *
 * Array based (T==cudaArray) textures instead use one of { tex1D, tex }support the full set of options:
 *   * cudaTextureAddressMode:
 *     - cudaTextureAddressMode::cudaAddressModeWrap     = 0
 *     - cudaTextureAddressMode::cudaAddressModeClamp    = 1
 *     - cudaTextureAddressMode::cudaAddressModeMirror   = 2
 *     - cudaTextureAddressMode::cudaAddressModeBorder   = 3
 *   * cudaTextureFilterMode:
 *     - cudaTextureFilterMode::cudaFilterModePoint      = 0
 *     - cudaTextureFilterMode::cudaFilterModeLinear     = 1
 *   * cudaTextureReadMode:
 *     - cudaTextureReadMode::cudaReadModeElementType    = 0
 *     - cudaTextureReadMode::cudaReadModeNormalizeFloat = 0
 *
 * Surfaces may only be bound to a cudaArray
 *
 *
 * */


template<typename T>
cudaResourceDesc _make_cudaResourceDesc(T *res, cudaResourceType restype=cudaResourceType::cudaResourceTypeLinear) {
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = restype;
    return resDesc;
}
template<> cudaResourceDesc _make_cudaResourceDesc<float>(float *res, cudaResourceType restype) {
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = restype;
    resDesc.res.linear.devPtr = res;
    resDesc.res.linear.desc.f = cudaChannelFormatKind::cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32; // bits per channel
    return resDesc;
}
template<> cudaResourceDesc _make_cudaResourceDesc<cudaArray>(cudaArray *res, cudaResourceType restype) {
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceType::cudaResourceTypeArray;
    resDesc.res.array.array = res;
    return resDesc;
}

void makeSurfObject(cudaSurfaceObject_t *surfobj, cudaArray *res) {
    // convenience function for initializing surface object
    // specify surface data
    struct cudaResourceDesc resDesc = _make_cudaResourceDesc<cudaArray>(res);

    // create
    checkCudaErrors(cudaCreateSurfaceObject(surfobj, &resDesc));
}

template<typename T>
void makeTexObject(cudaTextureObject_t *texobj, T *res, int dims,
        cudaTextureAddressMode addressmode, cudaTextureFilterMode filtermode,
        bool normalizedCoords=false, cudaTextureReadMode readmode=cudaTextureReadMode::cudaReadModeElementType,
        cudaResourceType restype=cudaResourceType::cudaResourceTypeLinear
) {
    // convenience function for initializing texture object

    // check arguments
    if (addressmode != cudaTextureAddressMode::cudaAddressModeClamp && addressmode != cudaTextureAddressMode::cudaAddressModeBorder ) {
        throw std::invalid_argument("invalid texture address mode given.");
    }
    if (filtermode != cudaTextureFilterMode::cudaFilterModePoint) { throw std::invalid_argument("invalid texture filter mode given."); }
    if (normalizedCoords) { throw std::invalid_argument("cannot use normalized coordinate addressing in linear textures."); }

    // specify texture data
    struct cudaResourceDesc resDesc = _make_cudaResourceDesc<T>(res, restype);
    resDesc.res.linear.sizeInBytes = dims*sizeof(T);

    // specify texture sampling params
    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    for (int i=0; i<dims; i++) {
        texDesc.addressMode[i] = addressmode;
    }
    texDesc.filterMode = filtermode;
    texDesc.readMode = readmode;
    texDesc.normalizedCoords = normalizedCoords;

    // create
    checkCudaErrors(cudaCreateTextureObject(texobj, &resDesc, &texDesc, 0));
}
template<> void makeTexObject<cudaArray>(cudaTextureObject_t *texobj, cudaArray *res, int dims,
        cudaTextureAddressMode addressmode, cudaTextureFilterMode filtermode,
        bool normalizedCoords, cudaTextureReadMode readmode, cudaResourceType restype
) {
    // convenience function for initializing texture object
    // specify texture data
    struct cudaResourceDesc resDesc = _make_cudaResourceDesc<cudaArray>(res);

    // specify texture sampling params
    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    for (int i=0; i<dims; i++) {
        texDesc.addressMode[i] = addressmode;
    }
    texDesc.filterMode = filtermode;
    texDesc.readMode = readmode;
    texDesc.normalizedCoords = normalizedCoords;

    // create
    checkCudaErrors(cudaCreateTextureObject(texobj, &resDesc, &texDesc, 0));

}

#endif // __MAKE_TEX_SURF_CUH__
