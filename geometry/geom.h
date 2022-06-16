#ifndef GEOM
#define GEOM

#include <array>
#include <vector>
#include <cuda_runtime.h>
#include "args.h"

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
    float SAD;
    float pixel_size;

    std::array<int, 2> fluence_map_dimension;
    std::array<int, 2> convolved_fluence_map_dimension;
    std::array<int, 2> extended_fluence_map_dimension;

    float* h_fluence_map;

    float* d_convolved_fluence_map;
    float* d_extended_fluence_map;

    beam();
    void convolve(FCBBkernel* kernel); // convolve from d_extended_fluence_map to d_convolved_fluence_map
};

void beams_init(std::vector<beam>& beams);
void test_convolve();
void host_convolve(float* h_convolved_fluence_map, \
    float* h_extended_fluence_map, float* convolution_kernel);

};

#endif