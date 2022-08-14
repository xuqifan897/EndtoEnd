#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <array>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

void fluence_map_init_water(beam& Beam, uint diam)
{
    // parameter diam specifies the diameter of the active fluence map
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = array<int, 2>( \
        {FM_dimension + 2 * FM_convolution_radius, FM_dimension + 2 * FM_convolution_radius});
    Beam.extended_fluence_map_dimension = array<int, 2>( \
        {FM_dimension + 4 * FM_convolution_radius, FM_dimension + 4 * FM_convolution_radius});
    
    size_t fluence_size = FM_dimension * FM_dimension;
    size_t convolved_size = Beam.convolved_fluence_map_dimension[0] * \
        Beam.convolved_fluence_map_dimension[1];
    size_t extended_size = Beam.extended_fluence_map_dimension[0] * \
        Beam.extended_fluence_map_dimension[1];
    
    float* h_extended_fluence_map = (float*)malloc(extended_size * sizeof(float));
    checkCudaErrors(cudaMalloc(&(Beam.d_convolved_fluence_map), convolved_size * sizeof(float)));
    checkCudaErrors(cudaMalloc(&(Beam.d_extended_fluence_map), extended_size * sizeof(float)));

    for (uint i=0; i<extended_size; i++)
        h_extended_fluence_map[i] = 0;

    uint margin = (Beam.extended_fluence_map_dimension[0] - diam) / 2;
    for (uint i=margin; i<margin + diam; i++)
    {
        uint I = i * Beam.extended_fluence_map_dimension[1];
        for (uint j=margin; j<margin + diam; j++)
        {
            uint idx = I + j;
            h_extended_fluence_map[idx] = 1;
        }
    }

    checkCudaErrors(cudaMemcpy(Beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_size*sizeof(float), cudaMemcpyHostToDevice));
    free(h_extended_fluence_map);
}

void E2E::test_FCBB_water_phantom(phantom& Phtm)
{
    // // phantom initialization
    // phantom Phtm;
    // phantom_init_default(Phtm);
    // Phtm.to_device();
    // Phtm.textureInit();

    // beam initialization
    beam Beam;
    Beam.zenith = PI / 2;
    Beam.azimuth = 0;
    Beam.SAD = get_args<float>("SAD") / 10;
    Beam.pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
    Beam.isocenter = Phtm.isocenter;
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = array<int, 2>({FM_dimension + 2 * FM_convolution_radius, \
        FM_dimension + 2 * FM_convolution_radius});
    Beam.extended_fluence_map_dimension = array<int, 2>({FM_dimension + 4 * FM_convolution_radius, \
        FM_dimension + 4 * FM_convolution_radius});
    Beam.FCBBinit(Phtm);

    // convolved_fluence_map, extended_fluence_map initialization
    fluence_map_init_water(Beam, 128);
    (*FCBB6MeV).d_conv_kernel_init();
    (*FCBB6MeV).texInit();
    Beam.convolve(FCBB6MeV);

    Beam.BEV_dose_forward(Phtm);
    Beam.PVCS_dose_forward(Phtm);

    uint volume = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_FCBB_PVCS_dose = (float*)malloc(volume * sizeof(float));
    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose, Beam.d_FCBB_PVCS_dose, \
        volume*sizeof(float), cudaMemcpyDeviceToHost));
    // string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/water_out/waterDose.dat"};
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/water_out/patientDose.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_FCBB_PVCS_dose, volume*sizeof(float));
    outFile.close();
}