#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <math.h>

#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;
namespace fs = boost::filesystem;


void phantom_log(phantom& Phtm)
{
    const array<int, 3>& phantom_dimension = Phtm.dimension;
    cout << "phantom shape: (" << phantom_dimension[0] << ", " << \
        phantom_dimension[1] << ", " << phantom_dimension[2] << ")" << endl;
    cout << "after pitch padding, the new shape is: (" << phantom_dimension[0] << \
        ", " << phantom_dimension[1] << ", " << Phtm.pitch << ")" << endl;
    cout << "The padding is to make the last dimension divisible by module " << \
        Phtm.pitch_module << ", for more efficient memory access." << endl;
    cout << "All outputs are in the new shape." << endl;
    cout << "isocenter: (" << Phtm.isocenter[0] << ", " << Phtm.isocenter[1] << ", " << Phtm.isocenter[2] << ")" << endl;
}


void kernel_log(FCBBkernel* kernel)
{
    cout << "The depth-dose kernel is a lookup table, with its depth from " \
        << (*kernel).min_depth << " to " << (*kernel).max_depth << " (cm), " \
        << (*kernel).num_depths << " points in total" << endl;
    cout << "Its parameters: A = " << (*kernel).A << ", B = " << (*kernel).B \
        << ", a = " << (*kernel).a << ", b = " << (*kernel).b << endl;
}

void beam_init(beam& Beam, phantom& Phtm)
{
    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    vector<float> fluence_map_dimension = get_args<vector<float>>("fluence-map-dimension");
    // convert mm to cm
    isocenter[0] /= 10;
    isocenter[1] /= 10;
    isocenter[2] /= 10;

    // randomly initialize zenith and azimuth angles
    uint mod = 10086;
    Beam.zenith = (float)(rand() % mod) / (mod - 1) * PI;
    Beam.azimuth = ((float)(rand() % mod) / (mod - 1) - 0.5) * 2 * PI;
    Beam.SAD = get_args<float>("SAD") / 10; // convert mm to cm
    Beam.pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = array<int, 2>({ \
        FM_dimension + 2 * FM_convolution_radius, FM_dimension + 2 * FM_convolution_radius});
    Beam.extended_fluence_map_dimension = array<int, 2>({ \
        FM_dimension + 4 * FM_convolution_radius, FM_dimension + 4 * FM_convolution_radius});
    Beam.isocenter = array<float, 3>({isocenter[0], isocenter[1], isocenter[2]});

    // log beam parameters
    cout << "zenith: " << Beam.zenith << " rad" << endl;
    cout << "azimuth: " << Beam.azimuth << " rad" << endl;
    cout << "SAD: " << Beam.SAD << " cm" << endl;
    cout << "pixel size: " << Beam.pixel_size << " cm" << endl;
    cout << "fluence map dimension: (" << Beam.fluence_map_dimension[0] << \
        ", " << Beam.fluence_map_dimension[1] << ")" << endl;
    cout << "convolved fluence map dimension: (" << Beam.convolved_fluence_map_dimension[0] << \
        ", " << Beam.convolved_fluence_map_dimension[1] << ")" << endl;
    cout << "extended fluence map dimension: (" << Beam.extended_fluence_map_dimension[0] << \
        ", " << Beam.extended_fluence_map_dimension[1] << ")" << endl;
    cout << "isocenter: (" << Beam.isocenter[0] << ", " << Beam.isocenter[1] << ", " << \
        Beam.isocenter[2] << ")" << endl;
    cout << "sampling range: (" << Beam.sampling_range[0] << ", " \
        << Beam.sampling_range[1] << ") cm" << endl;
    cout << "number of sampling points: " << Beam.sampling_points << endl;

    // cuda memory allocation
    auto& convolved_fluence_map_dimension = Beam.convolved_fluence_map_dimension;
    auto& extended_fluence_map_dimension = Beam.extended_fluence_map_dimension;
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map)), \
        convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_extended_fluence_map)), \
        extended_fluence_map_dimension[0]*extended_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map_grad)), \
        convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_fluence_grad)), \
        fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_element_wise_fluence_smoothness_loss)), \
        fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));

    // initialize extended fluence map randomly
    uint extended_dimension = FM_dimension + 4 * FM_convolution_radius;
    uint extended_fluence_map_size = extended_dimension * extended_dimension;
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size * sizeof(float));
    for (int i=0; i<extended_fluence_map_size; i++)
        h_extended_fluence_map[i] = 0;
    for (int i=0; i<FM_dimension; i++)
    {
        uint row_idx = (i + 2 * FM_convolution_radius) * extended_dimension;
        for (int j=0; j<FM_dimension; j++)
        {
            uint idx = row_idx + j + 2 * FM_convolution_radius;
            h_extended_fluence_map[idx] = (float)(rand() % mod) / (mod - 1);
        }
    }
    checkCudaErrors(cudaMemcpy(Beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
    free(h_extended_fluence_map);

    Beam.FCBBinit(Phtm);
}

float benchmark_conv(beam& Beam, FCBBkernel* kernel, uint rounds)
{
    float time = 0; // unit: millisecond.
    float milliseconds;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    for (int i=0; i<rounds; i++)
    {
        cudaEventRecord(start);
        Beam.convolve(kernel);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);

        cudaEventElapsedTime(&milliseconds, start, stop);
        cout << "conv benchmark round " << i << " time elapsed: " << milliseconds << endl;

        // skip the first round
        if (i==0)
            continue;
        time += milliseconds;
    }
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    return time / (rounds - 1);
}

float benchmark_convT(beam& Beam, FCBBkernel* kernel, uint rounds)
{
    // this function measures the back convolution time elapse
    float time = 0; // unit: milliseconds.
    float milliseconds; // store the time elapse for a single round
    uint mod = 10086; // module of random number

    // initialize convolved_fluence_map_grad randomly
    uint convolve_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolve_size = convolve_dimension * convolve_dimension;
    float* h_convolved_fluence_map_grad = (float*)malloc(convolve_size * sizeof(float));
    for (int i=0; i<convolve_size; i++)
        h_convolved_fluence_map_grad[i] = 0;
    for (int i=0; i<FM_dimension; i++)
    {
        uint row_idx = (i + FM_convolution_radius) * convolve_dimension;
        for (int j=0; j<FM_dimension; j++)
        {
            uint idx = row_idx + j + FM_convolution_radius;
            h_convolved_fluence_map_grad[idx] = (float)(rand() % mod) / (mod - 1);
        }
    }
    checkCudaErrors(cudaMemcpy(Beam.d_convolved_fluence_map_grad, \
        h_convolved_fluence_map_grad, convolve_size*sizeof(float), cudaMemcpyHostToDevice));
    free(h_convolved_fluence_map_grad);

    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    for (int i=0; i<rounds; i++)
    {
        cudaEventRecord(start);
        Beam.convolveT(kernel);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);

        cudaEventElapsedTime(&milliseconds, start, stop);
        cout << "convT benchmark round " << i << " time elapsed: " << milliseconds << endl;

        // skip the first round
        if (i == 0)
            continue;
        time += milliseconds;
    }
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    return time / (rounds - 1);
}

void benchmark_forward(beam& Beam, phantom& Phtm, FCBBkernel* kernel, \
    uint rounds, float* BEV_forward, float* PVCS_forward)
{
    float milliseconds_BEV;
    float milliseconds_PVCS;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));

    // initialize convolved_fluence_map
    uint convolve_dimension = FM_dimension + 2 * FM_convolution_radius;
    uint convolve_size = convolve_dimension * convolve_dimension;
    float* h_convolved_fluence_map = (float*)malloc(convolve_size * sizeof(float));
    for (int i=0; i<convolve_size; i++)
        h_convolved_fluence_map[i] = 0;
    uint mod = 10086;
    for (int i=0; i<FM_dimension; i++)
    {
        int row_idx = (i + FM_convolution_radius) * convolve_dimension;
        for (int j=0; j<FM_dimension; j++)
        {
            int idx = row_idx + j + FM_convolution_radius;
            h_convolved_fluence_map[idx] = (float)(rand() % mod) / (mod - 1);
        }
    }
    checkCudaErrors(cudaMemcpy(Beam.d_convolved_fluence_map, h_convolved_fluence_map, \
        convolve_size*sizeof(float), cudaMemcpyHostToDevice));
    free(h_convolved_fluence_map);

    *BEV_forward = 0;
    *PVCS_forward = 0;
    for (int i=0; i<rounds; i++)
    {   
        // set beamn angles
        Beam.zenith = (float)(rand() % mod) / (mod - 1) * PI;
        Beam.azimuth = ((float)(rand() % mod) / (mod - 1) - 0.5) * 2 * PI;

        cudaEventRecord(start);
        Beam.BEV_dose_forward(Phtm, kernel);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds_BEV, start, stop);

        cudaEventRecord(start);
        Beam.PVCS_dose_forward(Phtm);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds_PVCS, start, stop);

        cout << "round " << i << " BEV forward time: " << milliseconds_BEV \
            << ", PVCS forward time: " << milliseconds_PVCS << endl;;
        
        // skip the first round
        if (i == 0)
            continue;
        (*BEV_forward) += milliseconds_BEV;
        (*PVCS_forward) += milliseconds_PVCS;
    }
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    (*BEV_forward) /= (rounds - 1);
    (*PVCS_forward) /= (rounds - 1);
}

extern "C"
void writeFCBBPVCSDoseGradSurface(cudaSurfaceObject_t surface, float* input, \
    uint dim_x, uint dim_y, uint dim_z, cudaStream_t stream=0);

void benchmark_backward(beam& Beam, phantom& Phtm, FCBBkernel* kernel, \
    uint rounds, float* BEV_backward, float* PVCS_backward)
{
    float milliseconds_BEV;
    float milliseconds_PVCS;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    (*BEV_backward) = 0;
    (*PVCS_backward) = 0;

    // initialize FCBBPVCSDosegradSurface
    array<uint, 3> Phantom_shape {Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    uint Phantom_size {Phantom_shape[0] * Phantom_shape[1] * Phantom_shape[2]};
    float* h_FCBBPVCSDosegrad = (float*)malloc(Phantom_size * sizeof(float));
    uint mod = 10086;
    for (int i=0; i<Phantom_size; i++)
        h_FCBBPVCSDosegrad[i] = (float)(rand() % mod) / (mod - 1);
    float* d_FCBBPVCSDosegrad;
    checkCudaErrors(cudaMalloc((void**)(&d_FCBBPVCSDosegrad), Phantom_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_FCBBPVCSDosegrad, h_FCBBPVCSDosegrad, \
        Phantom_size*sizeof(float), cudaMemcpyHostToDevice));
    free(h_FCBBPVCSDosegrad);

    writeFCBBPVCSDoseGradSurface(beam::FCBB_PVCS_dose_grad_surface, \
        d_FCBBPVCSDosegrad, Phantom_shape[0], Phantom_shape[1], Phantom_shape[2]);

    for (int i=0; i<rounds; i++)
    {
        cudaEventRecord(start);
        Beam.PVCS_dose_backward(Phtm);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds_PVCS, start, stop);

        cudaEventRecord(start);
        Beam.BEV_dose_backward(Phtm, kernel);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds_BEV, start, stop);

        cout << "round " << i << " PVCS backward time: " << milliseconds_PVCS \
            << ", BEV backward time: " << milliseconds_BEV << endl;
        
        // skip the first round
        if (i == 0)
            continue;
        (*BEV_backward) += milliseconds_BEV;
        (*PVCS_backward) += milliseconds_PVCS;
    }

    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    (*BEV_backward) /= (rounds - 1);
    (*PVCS_backward) /= (rounds - 1);
}

int main(int argc, char** argv)
{
    // this function bechmarks the time for convolution step, forward step, and backward step
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }

    // phantom initialization
    cout << "\n\n\nPhantom initialization" << endl;
    phantom Phtm;
    phantom_init_default(Phtm);
    Phtm.to_device();
    Phtm.textureInit();
    phantom_log(Phtm);

    // kernel initialization
    cout << "\n\n\nKernel initialization" << endl;
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();
    (*kernel).texInit();
    kernel_log(kernel);

    // beam initialization
    // static initialization
    beam::FCBBStaticInit(Phtm);
    beam Beam;
    beam_init(Beam, Phtm);

    uint rounds = 11;
    float conv_bcmk_time = benchmark_conv(Beam, kernel, rounds);
    float convT_bcmk_time = benchmark_convT(Beam, kernel, rounds);
    float BEV_forward_time = 0;
    float PVCS_forward_time = 0;
    benchmark_forward(Beam, Phtm, kernel, rounds, &BEV_forward_time, &PVCS_forward_time);
    float BEV_backward_time = 0;
    float PVCS_backward_time = 0;
    benchmark_backward(Beam, Phtm, kernel, rounds, &BEV_backward_time, &PVCS_backward_time);

    cout << "\n\n" << endl;
    cout << "average convolution time: " << conv_bcmk_time << "ms." << endl;
    cout << "average transposed convolution time: " << convT_bcmk_time << "ms." << endl;
    cout << "average BEV forward time: " << BEV_forward_time << "ms." << endl;
    cout << "average PVCS forward time: " << PVCS_forward_time << "ms." << endl;
    cout << "average BEV backward time: " << BEV_backward_time << "ms." << endl;
    cout << "average PVCS backward time: " << PVCS_backward_time << "ms." << endl;
}