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
    // std::array<int, 3> dimension_org;
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
    float SAD; // in cm
    float pixel_size; // in cm

    std::array<int, 2> fluence_map_dimension;
    std::array<int, 2> convolved_fluence_map_dimension;
    std::array<int, 2> extended_fluence_map_dimension;

    float* h_fluence_map;

    float* d_convolved_fluence_map;
    float* d_extended_fluence_map;
    float* d_convolved_fluence_map_grad;
    float* d_fluence_grad;

    beam();
    // convolve from d_extended_fluence_map to d_convolved_fluence_map
    void convolve(FCBBkernel* kernel, cudaStream_t stream=0);

    /* compute the chain rule for gradient calculation. Essentially, 
    convolve with a convolution kernel rotated 180 degree. Here we 
    assume the convolution kernel is 180-degree rotation invariant */
    void convolveT(FCBBkernel* kernel, cudaStream_t stream=0);

    uint pitch_module;
    uint pitch;
    std::array<int, 3> dimension;
    std::array<float, 2> sampling_range;
    uint sampling_points;

    // logical order: (z, x, y), cudaExtent order: (y, x, z)
    cudaArray* FCBB_BEV_dose_array;
    cudaSurfaceObject_t FCBB_BEV_dose_surface;
    cudaTextureObject_t FCBB_BEV_dose_texture;

    // logical order: (x, y, z), cudaExtent order: (z, y, x)
    static bool FCBB_PVCS_dose_grad_init;
    static cudaArray* FCBB_PVCS_dose_grad_array;
    static cudaSurfaceObject_t FCBB_PVCS_dose_grad_surface;
    static cudaTextureObject_t FCBB_PVCS_dose_grad_texture;

    float* d_FCBB_PVCS_dose;

    void FCBBinit(phantom& Phtm);
    static void FCBBStaticInit(phantom& Phtm);

    void BEV_dose_forward(phantom& Phtm, FCBBkernel* FCBBkernel=FCBB6MeV, cudaStream_t stream=0);
    void PVCS_dose_forward(phantom& Phtm, cudaStream_t stream=0);
    static void calc_FCBB_PVCS_dose_grad(phantom& Phtm, float** d_elementWiseLoss, \
        float* d_PVCS_total_dose, cudaStream_t stream=0);

    /* Following FCBB_BEV_dose_array, d_FCBB_BEV_dose_grad 
    follows logical order: (z, x, y) */
    float* d_FCBB_BEV_dose_grad;
    void PVCS_dose_backward(phantom& Phtm, cudaStream_t stream=0);
    void BEV_dose_backward(phantom& Phtm, FCBBkernel* kernel=FCBB6MeV, cudaStream_t stream=0);
    void fluence_map_update(uint idx, float* d_norm_final, float* d_squared_grad, float step_size, cudaStream_t stream=0);
};

void beams_init(std::vector<beam>& beams);
void test_convolve();
void test_convolveT();
void host_convolve(float* h_convolved_fluence_map, \
    float* h_extended_fluence_map, float* convolution_kernel, \
    uint convolved_fluence_map_size, \
    uint extended_fluence_map_size, uint source_prepend=0);

void test_volume_rendering();
void test_BEV_dose_forward();
void test_PVCS_surface();
void test_PVCS_dose_forward();
void test_FCBB_water_phantom();
void test_calc_FCBB_PVCS_dose_grad(std::vector<beam>& beams, phantom& Phtm);
void test_FCBB_PVCS_backward(std::vector<beam>& beams, phantom& Phtm);
void test_minus_coordinates_of_texture_memory_out_of_curiosity();
void test_FCBB_BEV_backward(std::vector<beam>& beams, phantom& Phtm);
void test_fluence_map_update(std::vector<beam>& beams);

// for debug purposes
extern float* HU_debug;
extern float* dose_debug;
extern bool* valid_debug;
};

#endif