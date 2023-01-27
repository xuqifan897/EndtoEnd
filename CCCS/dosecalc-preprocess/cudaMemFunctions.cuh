#ifndef __CUDAMEMFUNCTIONS_H__
#define __CUDAMEMFUNCTIONS_H__

#include <cuda_runtime.h>
#include "helper_cuda.h"
#include "helper_math.h"
#include "DoseCalcIO/volume.h"
#include "DoseCalcIO/roi.h"
#include "DoseCalcIO/ctlut.h"

int read_dicom(FloatVolume& dens, FloatVolume& ctdata, float iso_voxel, CTLUT* ctlut=nullptr, bool verbose=false);
void cudaHU2dens(FloatVolume&, ShortVolume&, float, bool );
void cudaCreateTexIso( FloatVolume& dens, FloatVolume& data, float iso_size, CTLUT* ctlut=nullptr, bool hu2dens=true);
Volume<uint8_t> generateContourMask(StructureSet&, FrameOfReference, FrameOfReference, void* texRay=NULL);
void findFluenceProjection( std::vector<float>& fluence,
                            const Volume<uint8_t>& ptv_mask,
                            float3 isocenter,
                            float3 source,
                            uint2 fmap_size,
                            float2 beamlet_size,
                            float azimuth, float zenith, float coll, int verbose);

#endif // __CUDAMEMFUNCTIONS_H__
