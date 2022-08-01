#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

void E2E::module_test_convolve(beam& Beam)
{
    uint convolved_fluence_map_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint kernel_dimension = 2 * FM_convolution_radius;
    uint extended_fluence_map_dimension = FM_dimension + 4 * FM_convolution_radius;

    // initialize extended_fluence_map
    uint extended_fluence_map_size = extended_fluence_map_dimension * \
        extended_fluence_map_dimension;
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size*sizeof(float));
    uint norm = 65536;
    for (uint i=0; i<extended_fluence_map_size; i++)
        h_extended_fluence_map[i] = (float)(rand() % norm) / (norm - 1);
    
    checkCudaErrors(cudaMemcpy(Beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    
    // initialize convolution kernel
    uint convolution_kernel_size = kernel_dimension * kernel_dimension;
    float* h_convolution_kernel = (float*)malloc(convolution_kernel_size*sizeof(float));
    for (uint i=0; i<convolution_kernel_size; i++)
        h_convolution_kernel[i] = (float)(rand() % norm) / (norm - 1);

    checkCudaErrors(cudaMemcpy(FCBB6MeV->d_convolution_kernel, h_convolution_kernel, \
        convolution_kernel_size*sizeof(float), cudaMemcpyHostToDevice));
    
    // initialize output
    uint convolved_fluence_map_size = convolved_fluence_map_dimension * \
        convolved_fluence_map_dimension;
    float* h_convolved_fluence_map = (float*)malloc(convolved_fluence_map_size*sizeof(float));
    float* h_convolved_fluence_map_d = (float*)malloc(convolved_fluence_map_size*sizeof(float));

    // device compute
    Beam.convolve(FCBB6MeV, 0);
    checkCudaErrors(cudaMemcpy(h_convolved_fluence_map_d, Beam.d_convolved_fluence_map, \
        convolved_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
    
    // host compute
    constexpr uint source_prepend = 1;
    for (uint i=0; i<convolved_fluence_map_dimension; i++)
    {
        for (uint j=0; j<convolved_fluence_map_dimension; j++)
        {
            uint convolved_idx = i * convolved_fluence_map_dimension + j;
            h_convolved_fluence_map[convolved_idx] = 0;
            for (uint k=0; k<kernel_dimension; k++)
            {
                for (uint l=0; l<kernel_dimension; l++)
                {
                    uint kernel_idx = k * kernel_dimension + l;
                    uint extended_idx = (i + k + source_prepend) * extended_fluence_map_dimension + j + l + source_prepend;
                    h_convolved_fluence_map[convolved_idx] += h_convolution_kernel[kernel_idx] * h_extended_fluence_map[extended_idx];
                }
            }
        }
    }

    // compare
    float mse = 0;
    for (uint i=0; i<convolved_fluence_map_size; i++)
    {
        float diff = h_convolved_fluence_map_d[i] - h_convolved_fluence_map[i];
        mse += diff * diff;
    }
    mse = sqrt(mse / convolved_fluence_map_size);

    float base = 0;
    for (uint i=0; i<convolved_fluence_map_size; i++)
        base += h_convolved_fluence_map[i] * h_convolved_fluence_map[i];
    base = sqrt(base / convolved_fluence_map_size);

    cout << "The MSE value between host and device results is " << mse << ", while the scale of host result is " << base << endl;

    // cleanup
    free(h_extended_fluence_map);
    free(h_convolution_kernel);
    free(h_convolved_fluence_map);
    free(h_convolved_fluence_map_d);
}

void test_modules_beam_init(beam& Beam, phantom& Phtm)
{
    Beam.zenith = PI / 2;
    Beam.azimuth = 0;
    Beam.SAD = get_args<float>("SAD") / 10;
    Beam.pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.isocenter = Phtm.isocenter;
    Beam.convolved_fluence_map_dimension = array<int, 2>({FM_dimension + 2 * FM_convolution_radius, \
        FM_dimension + 2 * FM_convolution_radius});
    Beam.extended_fluence_map_dimension = array<int, 2>({FM_dimension + 4 * FM_convolution_radius, \
        FM_dimension + 4 * FM_convolution_radius});
    Beam.FCBBinit(Phtm);
    beam::FCBBStaticInit(Phtm);

    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map)), \
        Beam.convolved_fluence_map_dimension[0]*Beam.convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_extended_fluence_map)), \
        Beam.extended_fluence_map_dimension[0]*Beam.extended_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map_grad)), \
        Beam.convolved_fluence_map_dimension[0]*Beam.convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_fluence_grad)), \
        Beam.fluence_map_dimension[0]*Beam.fluence_map_dimension[1]*sizeof(float)));
}

void E2E::test_modules(phantom& Phtm)
{
    // beam initialization
    beam Beam;
    test_modules_beam_init(Beam, Phtm);

    srand(1008611);
    // module_test_convolve(Beam);
    module_test_BEV_dose_forward(Beam, Phtm);
        // module_test_host_triliner(Phtm);
        // module_test_host_linear();
}