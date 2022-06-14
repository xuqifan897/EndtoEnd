#ifndef GEOM
#define GEOM

#include <array>
#include <cuda_runtime.h>

namespace E2E
{

class phantom
{
public:
    std::array<int, 3> dimension;
    std::array<float, 3> isocenter; // in cm
    float voxelSize; // in cm, isotropic phantom is assumed

    int pitch_module;
    int pitch; // pitch in number of floats
    bool pitchPadding; // whether the arrays are pitched;

    float* h_HU; // the water HU value is normalized to 1
    float* h_PTVweight;
    float* h_PTVtarget;
    float* h_OARweight;
    float* h_OARtarget;
    float* h_Dose;

    cudaArray* d_HU;
    float* d_PTVweight;
    float* d_PTVtarget;
    float* d_OARweight;
    float* d_OARtarget;
    float* d_Dose;

    cudaTextureObject_t tex;
    bool texInit;

    phantom();
    ~phantom();
    void pitchPad();
    void Dose_init();
    void to_device();
    void DoseToHost();
    void textureInit();
    void textureDecon();
};

int phantom_init_default(phantom& Phtm);
void runTest(phantom& Phtm);

class beam
{
public:
    float zenith;
    float azimuth;
    std::array<float, 2> fluence_size;
    std::array<float, 2> padded_fluence_size;

    float* h_fluence_map;
    float* h_padded_fluence_map;

    float* d_fluence_map;
    float* d_padded_fluence_map;
};

};

#endif