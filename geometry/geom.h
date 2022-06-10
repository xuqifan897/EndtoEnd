#ifndef GEOM
#define GEOM

#include <array>

namespace E2E
{

class phantom
{
public:
    std::array<int, 3> dimension;
    std::array<float, 3> isocenter; // in cm
    float voxelSize; // in cm, isotropic phantom is assumed
    float* h_HU; // the water HU value is normalized to 1
    float* h_PTVweight;
    float* h_PTVtarget;
    float* h_OARweight;
    float* h_OARtarget;

    phantom();
    ~phantom();
    phantom(phantom& old);
    phantom(phantom&& old);
};

int phantom_init_default(phantom& Phtm);

};

#endif