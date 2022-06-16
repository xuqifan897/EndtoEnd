#include <iostream>
#include <fstream>
#include <vector>
#include <array>
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

    this->h_fluence_map = nullptr;
    this->d_convolved_fluence_map = nullptr;
    this->d_extended_fluence_map = nullptr;
}

void E2E::beams_init(vector<beam>& beams)
{
    float SAD = get_args<float>("SAD") / 10; // mm to cm
    float fluence_map_pixel_size = get_args<float>("fluence-map-pixel-size");
    vector<int> fluence_map_convolution_radius = \
        get_args<vector<int>>("fluence-map-convolution-radius");
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
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_convolved_fluence_map)), \
            convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_extended_fluence_map)), \
            extended_fluence_map_dimension[0]*extended_fluence_map_dimension[1]*sizeof(float)));
    }
}


extern "C"
void convolve_kernel(dim3 gridSize, dim3 blockSize, float* convolvedFluenceMap, \
    float* extendedFluenceMap, float* convolutionKernel, int ConvRad, int globalPitch);

void beam::convolve(FCBBkernel* kernel)
{
    // float* d_convolution_kernel = (*kernel).d_convolution_kernel;
    int blockS = 16;
    dim3 blockSize(blockS, blockS);
    dim3 gridSize(this->convolved_fluence_map_dimension[0]/blockS, \
        this->convolved_fluence_map_dimension[1]/blockS);
    int globalPitch = FM_dimension + 4 * FM_convolution_radius;
    convolve_kernel(gridSize, blockSize, this->d_convolved_fluence_map, this->d_extended_fluence_map, \
        (*kernel).d_convolution_kernel, FM_convolution_radius, globalPitch);
}


void E2E::test_convolve()
{
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();

    // // for debug purposes
    // uint kernel_size = (*kernel).kernel_pitch;
    // float* h_convolution_kernel = (float*)malloc(kernel_size*kernel_size*sizeof(float));
    // checkCudaErrors(cudaMemcpy(h_convolution_kernel, (*kernel).d_convolution_kernel, \
    //     kernel_size*kernel_size*sizeof(float), cudaMemcpyDeviceToHost));
    // string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/convolutionKernel1.dat"};
    // ofstream outFile(outputPath);
    // outFile.write((char*)h_convolution_kernel, kernel_size*kernel_size*sizeof(float));
    // outFile.close();
    // exit;

    beam Beam;
    array<int, 2> convolved_fluence_map_dimension({FM_dimension+2*FM_convolution_radius, \
        FM_dimension+2*FM_convolution_radius});
    array<int, 2> extended_fluence_map_dimension({FM_dimension+4*FM_convolution_radius, \
        FM_dimension+4*FM_convolution_radius});
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
    host_convolve(h_convolved_fluence_map, h_extended_fluence_map, h_convolution_kernel);
    cout << "host code end" << endl;
    outputPath = "/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/convTest_host.dat";
    outFile.open(outputPath);
    outFile.write((char*)h_convolved_fluence_map, convolved_size*sizeof(float));
    outFile.close();
}


void E2E::host_convolve(float* h_convolved_fluence_map, \
    float* h_extended_fluence_map, float* convolution_kernel)
{
    uint kernelSize = 2 * FM_convolution_radius;
    uint convolved_fluence_map_size = FM_dimension + 2 * FM_convolution_radius;
    uint extended_fluence_map_size = FM_dimension + 4 * FM_convolution_radius;
    for (uint i=0; i<convolved_fluence_map_size; i++)
    {
        for (uint j=0; j<convolved_fluence_map_size; j++)
        {   
            size_t convolved_fluence_map_idx = i * convolved_fluence_map_size + j;
            h_convolved_fluence_map[convolved_fluence_map_idx] = 0;
            for (uint k=0; k<kernelSize-1; k++)
            {
                for (uint l=0; l<kernelSize-1; l++)
                {
                    size_t extended_fluence_map_idx = (i + k + 1) * extended_fluence_map_size + j + l + 1;
                    size_t convolution_kernel_idx = k * kernelSize + l;
                    h_convolved_fluence_map[convolved_fluence_map_idx] += \
                        h_extended_fluence_map[extended_fluence_map_idx] * convolution_kernel[convolution_kernel_idx];
                }
            }
        }
    }
}