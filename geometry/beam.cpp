#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <math.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"

using namespace E2E;
using namespace std;

beam::beam()
{
    this->zenith = 0;
    this->azimuth = 0;
    this->SAD = 0;
    this->pixel_size = 0;

    this->fluence_map_dimension = array<int, 2>({0, 0});
    this->convolved_fluence_map_dimension = array<int, 2>({0, 0});
    this->extended_fluence_map_dimension = array<int, 2>({0, 0});
    this->isocenter = array<float, 3>({0, 0});

    this->h_fluence_map = nullptr;
    this->d_convolved_fluence_map = nullptr;
    this->d_extended_fluence_map = nullptr;
    this->d_convolved_fluence_map_grad = nullptr;
    this->d_fluence_grad = nullptr;

    this->pitch_module = E2E::pitch_module;
    vector<int> dimension_ = get_args<vector<int>>("phantom-dimension");
    this->dimension = array<int, 3>({dimension_[0], dimension_[1], dimension_[2]});
    this->pitch = ceil((float)(this->dimension[2]) / this->pitch_module) * this->pitch_module;
    vector<float> sampling_range_ = get_args<vector<float>>("fluence-map-sampling-range");
    this->sampling_range = array<float, 2>({sampling_range_[0] / 10, sampling_range_[1] / 10});
    this->sampling_points = get_args<int>("fluence-map-sampling-points");

    this->FCBB_BEV_dose_array = 0;
    this->FCBB_BEV_dose_surface = 0;
    this->FCBB_BEV_dose_texture = 0;

    // // The following member variables have been changed to static
    // this->FCBB_PVCS_dose_grad_array = 0;
    // this->FCBB_PVCS_dose_grad_surface = 0;
    // this->FCBB_PVCS_dose_grad_texture = 0;
    
    this->d_FCBB_PVCS_dose = nullptr;
    this->d_FCBB_BEV_dose_grad = nullptr;
}

void E2E::beams_init(vector<beam>& beams)
{
    float SAD = get_args<float>("SAD") / 10; // mm to cm
    float fluence_map_pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
    vector<int> fluence_map_convolution_radius = \
        get_args<vector<int>>("fluence-map-convolution-radius");
    
    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    for (uint i=0; i<isocenter.size(); i++)
        isocenter[i] /= 10; // convert mm to cm
    
    if (fluence_map_convolution_radius[0]!=E2E::FM_convolution_radius || \
        fluence_map_convolution_radius[1]!=E2E::FM_convolution_radius)
    {
        cout << "Sorry, we only support fluence map \
            convolution radius of " << E2E::FM_convolution_radius << " at this time" << endl;
        exit;
    }
    vector<int> fluence_map_dimension = get_args<vector<int>>("fluence-map-dimension");
    if (fluence_map_dimension[0]!=E2E::FM_dimension || fluence_map_dimension[1]!=E2E::FM_dimension)
    {
        cout << "Sorry, we only support fluence map \
            dimension of " << E2E::FM_dimension << " at this time" << endl;
        exit;
    }

    vector<int> convolved_fluence_map_dimension{FM_dimension+2*FM_convolution_radius, \
        FM_dimension+2*FM_convolution_radius};
    vector<int> extended_fluence_map_dimension{FM_dimension+4*FM_convolution_radius, \
        FM_dimension+4*FM_convolution_radius};

    string beam_angle_config_path = get_args<string>("beam-angle-config-path");
    ifstream inFile(beam_angle_config_path);
    if (!inFile.is_open())
    {
        cout << "Could not open this file: " << beam_angle_config_path << endl;
        exit;
    }
    vector<string> lines;
    string line;
    while(getline(inFile, line))
        lines.push_back(line);
    inFile.close();
    
    for (int i=0; i<lines.size(); i++)
    {
        line = lines[i];
        vector<float> line_parsed;
        string item;
        stringstream ss(line);
        while(getline(ss, item, ','))
            line_parsed.push_back(stof(item));
        beams.push_back(beam());
        auto& new_beam = beams.back();
        new_beam.zenith = line_parsed[0];
        new_beam.azimuth = line_parsed[1];
        new_beam.SAD = SAD;
        new_beam.pixel_size = fluence_map_pixel_size;
        new_beam.fluence_map_dimension = array<int, 2>({fluence_map_dimension[0], \
            fluence_map_dimension[1]});
        new_beam.convolved_fluence_map_dimension = array<int, 2>({convolved_fluence_map_dimension[0], \
            convolved_fluence_map_dimension[1]});
        new_beam.extended_fluence_map_dimension = array<int, 2>( \
            {extended_fluence_map_dimension[0], extended_fluence_map_dimension[1]});
        new_beam.isocenter = array<float, 3>({isocenter[0], isocenter[1], isocenter[2]});
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_convolved_fluence_map)), \
            convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_extended_fluence_map)), \
            extended_fluence_map_dimension[0]*extended_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_convolved_fluence_map_grad)), \
            convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_fluence_grad)), \
            fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));
    }
}


extern "C"
void convolve_kernel(dim3 gridSize, dim3 blockSize, cudaStream_t stream, float* convolvedFluenceMap, \
    float* extendedFluenceMap, float* convolutionKernel, int ConvRad, int globalPitch,  
    int target_prepend=0, int source_prepend=0, int kernel_prepend=0);

void beam::convolve(FCBBkernel* kernel, cudaStream_t stream)
{
    // float* d_convolution_kernel = (*kernel).d_convolution_kernel;
    int blockS = 16;
    dim3 blockSize(blockS, blockS);
    dim3 gridSize(this->convolved_fluence_map_dimension[0]/blockS, \
        this->convolved_fluence_map_dimension[1]/blockS);
    int globalPitch = FM_dimension + 4 * FM_convolution_radius;
    convolve_kernel(gridSize, blockSize, stream, this->d_convolved_fluence_map, this->d_extended_fluence_map, \
        (*kernel).d_convolution_kernel, FM_convolution_radius, globalPitch, 0, 1, 0);
}

void beam::convolveT(FCBBkernel* kernel, cudaStream_t stream)
{   
    if (this->d_convolved_fluence_map_grad == nullptr)
    {
        cout << "member d_convolved_fluence_map_grad is not initialized." << endl;
        exit;
    }
    if (this->d_fluence_grad == nullptr)
    {
        uint fluence_grad_size = FM_dimension * FM_dimension;
        checkCudaErrors(cudaMalloc((void**)&(this->d_fluence_grad), fluence_grad_size * sizeof(float)));
    }
    uint blockS = 16;
    dim3 blockSize(blockS, blockS);
    dim3 gridSize(FM_dimension/blockS, FM_dimension/blockS);
    uint d_convolved_fluence_map_pitch = FM_dimension + 2 * FM_convolution_radius;
    convolve_kernel(gridSize, blockSize, stream, this->d_fluence_grad, this->d_convolved_fluence_map_grad, \
        (*kernel).d_convolution_kernel, FM_convolution_radius, d_convolved_fluence_map_pitch, 0, 0, 0);
}

void E2E::test_convolve()
{
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();

    // output the convolution kernel
    uint kernelSize = 4 * FM_convolution_radius * FM_convolution_radius;
    float* h_convKernel = (float*)malloc(kernelSize * sizeof(float));
    checkCudaErrors(cudaMemcpy(h_convKernel, (*kernel).d_convolution_kernel, \
        kernelSize * sizeof(float), cudaMemcpyDeviceToHost));
    string FCBBconvkernelPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/convKernelEven.dat"};
    ofstream outFile_(FCBBconvkernelPath);
    outFile_.write((char*)h_convKernel, kernelSize * sizeof(float));
    outFile_.close();

    beam Beam;
    array<int, 2> convolved_fluence_map_dimension({FM_dimension+2*FM_convolution_radius, \
        FM_dimension+2*FM_convolution_radius});
    array<int, 2> extended_fluence_map_dimension({FM_dimension+4*FM_convolution_radius, \
        FM_dimension+4*FM_convolution_radius});
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = convolved_fluence_map_dimension;
    Beam.extended_fluence_map_dimension = extended_fluence_map_dimension;
    
    size_t convolved_size = convolved_fluence_map_dimension[0] * convolved_fluence_map_dimension[1];
    size_t extended_size = extended_fluence_map_dimension[0] * extended_fluence_map_dimension[1];
    size_t kernel_size = 4 * FM_convolution_radius * FM_convolution_radius;
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map)), convolved_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_extended_fluence_map)), extended_size * sizeof(float)));
    
    float* h_convolved_fluence_map = (float*)malloc(convolved_size * sizeof(float));
    float* h_extended_fluence_map = (float*)malloc(extended_size * sizeof(float));
    float* h_convolution_kernel = (float*)malloc(kernel_size*sizeof(float));
    
    string extendedPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/extend.dat"};
    ifstream inFile(extendedPath);
    inFile.read((char*)h_extended_fluence_map, extended_size*sizeof(float));
    inFile.close();

    string convolutionKernelPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/kernel.dat"};
    inFile.open(convolutionKernelPath);
    inFile.read((char*)h_convolution_kernel, kernel_size*sizeof(float));
    inFile.close();
    
    checkCudaErrors(cudaMemcpy(Beam.d_extended_fluence_map, h_extended_fluence_map, \
        extended_size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy((*kernel).d_convolution_kernel, h_convolution_kernel, \
        kernel_size*sizeof(float), cudaMemcpyHostToDevice));
    
    Beam.convolve(kernel);
    checkCudaErrors(cudaMemcpy(h_convolved_fluence_map, Beam.d_convolved_fluence_map, \
        convolved_size*sizeof(float), cudaMemcpyDeviceToHost));
    
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/convTest.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_convolved_fluence_map, convolved_size*sizeof(float));
    outFile.close();

    // for debug purposes
    cout << "host code begin" << endl;
    host_convolve(h_convolved_fluence_map, h_extended_fluence_map, h_convolution_kernel, \
        convolved_fluence_map_dimension[0], extended_fluence_map_dimension[0], 1);
    cout << "host code end" << endl;
    outputPath = "/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/convTest_host.dat";
    outFile.open(outputPath);
    outFile.write((char*)h_convolved_fluence_map, convolved_size*sizeof(float));
    outFile.close();
}

void E2E::test_convolveT()
{
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();

    beam Beam;
    array<int, 2> convolved_fluence_map_dimension({FM_dimension+2*FM_convolution_radius, \
        FM_dimension+2*FM_convolution_radius});
    array<int, 2> extended_fluence_map_dimension({FM_dimension+4*FM_convolution_radius, \
        FM_dimension+4*FM_convolution_radius});
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = convolved_fluence_map_dimension;
    Beam.extended_fluence_map_dimension = extended_fluence_map_dimension;

    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_convolved_fluence_map_grad)), \
        convolved_fluence_map_dimension[0] * convolved_fluence_map_dimension[1] * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&(Beam.d_fluence_grad)), \
        FM_dimension*FM_dimension*sizeof(float)));

    float* h_convolved_fluence_map_grad = (float*)malloc(convolved_fluence_map_dimension[0] * \
        convolved_fluence_map_dimension[1] * sizeof(float));
    float* h_fluence_map_grad = (float*)malloc(FM_dimension * FM_dimension * sizeof(float));
    float* h_convolution_kernel = (float*)malloc((*kernel).kernel_pitch * \
        (*kernel).kernel_pitch * sizeof(float));
    
    string convolved_FM_grad_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/"
        "convolved_FM_grad.dat"};
    string kernelTpath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/kernelT.dat"};
    ifstream inFile(convolved_FM_grad_path);
    inFile.read((char*)h_convolved_fluence_map_grad, convolved_fluence_map_dimension[0] * \
        convolved_fluence_map_dimension[1] * sizeof(float));
    inFile.close();
    inFile.open(kernelTpath);
    inFile.read((char*)h_convolution_kernel, (*kernel).kernel_pitch * \
        (*kernel).kernel_pitch * sizeof(float));
    inFile.close();

    checkCudaErrors(cudaMemcpy(Beam.d_convolved_fluence_map_grad, h_convolved_fluence_map_grad, \
        convolved_fluence_map_dimension[0] * \
        convolved_fluence_map_dimension[1] * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy((*kernel).d_convolution_kernel, h_convolution_kernel, \
        (*kernel).kernel_pitch * (*kernel).kernel_pitch * sizeof(float),
        cudaMemcpyHostToDevice));
    Beam.convolveT(kernel);
    checkCudaErrors(cudaMemcpy(h_fluence_map_grad, Beam.d_fluence_grad, \
        FM_dimension*FM_dimension*sizeof(float), cudaMemcpyDeviceToHost));
    string d_FM_grad_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/d_FM_grad.dat"};
    ofstream outFile(d_FM_grad_path);
    outFile.write((char*)h_fluence_map_grad, FM_dimension*FM_dimension*sizeof(float));
    outFile.close();

    host_convolve(h_fluence_map_grad, h_convolved_fluence_map_grad, h_convolution_kernel, \
        FM_dimension, convolved_fluence_map_dimension[0]);
    string h_FM_grad_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/h_FM_grad.dat"};
    outFile.open(h_FM_grad_path);
    outFile.write((char*)h_fluence_map_grad, FM_dimension*FM_dimension*sizeof(float));
}

void E2E::host_convolve(float* h_convolved_fluence_map, \
    float* h_extended_fluence_map, float* convolution_kernel, \
    uint convolved_fluence_map_size, uint extended_fluence_map_size, uint source_prepend)
{
    uint kernelSize = 2 * FM_convolution_radius;
    // uint convolved_fluence_map_size = FM_dimension + 2 * FM_convolution_radius;
    // uint extended_fluence_map_size = FM_dimension + 4 * FM_convolution_radius;
    for (uint i=0; i<convolved_fluence_map_size; i++)
    {
        for (uint j=0; j<convolved_fluence_map_size; j++)
        {   
            size_t convolved_fluence_map_idx = i * convolved_fluence_map_size + j;
            h_convolved_fluence_map[convolved_fluence_map_idx] = 0;
            for (uint k=0; k<kernelSize; k++)
            {
                for (uint l=0; l<kernelSize; l++)
                {
                    size_t extended_fluence_map_idx = (i + k + source_prepend) * \
                        extended_fluence_map_size + j + l + source_prepend;
                    size_t convolution_kernel_idx = k * kernelSize + l;
                    h_convolved_fluence_map[convolved_fluence_map_idx] += \
                        h_extended_fluence_map[extended_fluence_map_idx] * convolution_kernel[convolution_kernel_idx];
                }
            }
        }
    }
}