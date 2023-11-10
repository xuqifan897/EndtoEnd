#include "calc_header.h"
#include "argparse.h"

void old::constantsInit(CONSTANTS* host)
{
    if (host == nullptr)
        host = new CONSTANTS{};

    const auto& dicomVolumeDimension = dev::getarg<std::vector<int>>("dicomVolumeDimension");
    (host)->size = make_uint3(dicomVolumeDimension[0], dicomVolumeDimension[1], dicomVolumeDimension[2]);
    const auto& voxelSize = dev::getarg<std::vector<float>>("voxelSize");
    (host)->voxel = make_float3(voxelSize[0], voxelSize[1], voxelSize[2]);
    const auto& dicomVolumeStartCoords = dev::getarg("dicomVolumeStartCoords");
    (host)->start = make_float3(dicomVolumeStartCoords[0], dicomVolumeStartCoords[1], dicomVolumeStartCoords[2]);
    const auto& doseBoundingBoxStartIndices = dev::getarg<std::vector<int>>("doseBoundingBoxStartIndices");
    (host)->calc_bbox_start = make_uint3(doseBoundingBoxStartIndices[0], doseBoundingBoxStartIndices[1], doseBoundingBoxStartIndices[2]);
    const auto& doseBoundingBoxDimensions = dev::getarg<std::vector<int>>("doseBoundingBoxDimensions");
    (host)->calc_bbox_size = make_uint3(doseBoundingBoxDimensions[0], doseBoundingBoxDimensions[1], doseBoundingBoxDimensions[2]);
    const auto& REVConvolutionArrayDimensions = dev::getarg<std::vector<int>>("REVConvolutionArrayDimensions");
    (host)->max_rev_size = make_uint3(REVConvolutionArrayDimensions[0], REVConvolutionArrayDimensions[1], REVConvolutionArrayDimensions[2]);
    float convlat = dev::getarg<float>("convlat");
    (host)->rev_latspacing = convlat;
    (host)->rev_longspacing = dev::getarg<float>("convstep");
    (host)->kernel_extent = dev::getarg<float>("kernelExtent");
    (host)->ss_factor = 1;
    (host)->nphi = static_cast<uint>(dev::getarg<int>("nphi"));
    (host)->ntheta = static_cast<uint>(dev::getarg<int>("ntheta"));
    (host)->nradii = static_cast<uint>(dev::getarg<int>("nradii"));
    (host)->penumbra = dev::getarg<float>("penumbra");
    (host)->beam_count = dev::getarg<int>("beamCount");
    (host)->beam_spec = std::string("spec_6mv.spec");
}