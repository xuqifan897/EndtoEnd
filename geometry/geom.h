#ifndef GEOM
#define GEOM

#include <array>


namespace E2E
{
    class phantom
    {
    public:
        std::array<int, 3> dimension;
        float voxelSize; // in cm
        std::array<float, 3> isoCenter; // in cm
        float* d_HU;
        float* d_PTVtarget;
        float* d_PTVweight;
        // float* d_PTVtarget;
        float* d_OARweight;
        float* d_OARtarget;

        phantom();
        ~phantom();
        phantom(phantom& old);
        phantom(phantom&& old);
    };

    int phantom_init_default(phantom& target);
};

#endif