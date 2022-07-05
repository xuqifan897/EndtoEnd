#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>
#include <math.h>
#include "geom.h"
#include "args.h"

using namespace E2E;
using namespace std;

void beam::FCBBinit(phantom& Phtm)
{
    // const cudaExtent volumeSize_ = make_cudaExtent( \
    //     this->pitch, this->dimension[1], this->dimension[0]);

    /* here, due to the radiological path calculation is along z direction, 
    we want x and y to be contiguous to facilitate data read and write 
    logical: (z, x, y), for cudaExtent: (y, x, z)*/
    const cudaExtent volumeSize_ = make_cudaExtent(this->convolved_fluence_map_dimension[1], \
        this->convolved_fluence_map_dimension[0], this->sampling_points);
    
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    checkCudaErrors(cudaMalloc3DArray(&(this->FCBB_BEV_dose_array), \
        &channelDesc, volumeSize_, cudaArraySurfaceLoadStore));

    cudaResourceDesc surfRes;
    memset(&surfRes, 0, sizeof(cudaResourceDesc));
    surfRes.resType = cudaResourceTypeArray;
    surfRes.res.array.array = this->FCBB_BEV_dose_array;
    checkCudaErrors(cudaCreateSurfaceObject(&(this->FCBB_BEV_dose_surface), &surfRes));

    cudaResourceDesc texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array = this->FCBB_BEV_dose_array;

    cudaTextureDesc texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));

    texDescr.normalizedCoords = true; // access with normalized texture coordinates
    texDescr.filterMode = cudaFilterModeLinear; // linear interpolation
    // wrap texture coordinates
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.addressMode[1] = cudaAddressModeBorder;
    texDescr.addressMode[2] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;

    checkCudaErrors(cudaCreateTextureObject(&(this->FCBB_BEV_dose_texture), \
        &texRes, &texDescr, NULL));

    checkCudaErrors(cudaMalloc(&d_FCBB_PVCS_dose, \
    Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch * sizeof(float)));
    
    /* Ensure that the phantom voxel size is larger than the largest 
    dimension of a beam cell, so that during gradient calculation, each beam 
    point can have at most 8 sources*/
    
    float pixel_size_at_the_furthest_range = this->pixel_size * this->sampling_range[1] / this->SAD;
    float vector_x[3]{pixel_size_at_the_furthest_range, 0, 0};
    float vector_y[3]{0, pixel_size_at_the_furthest_range, 0};

    // float vector_z_prime[3]{-(this->dimension[0] / 2 - 1) * this->pixel_size, \
    //     -(this->dimension[1] / 2 - 1) * this->pixel_size, -this->SAD};
    float vector_z_prime[3]{-(float)(this->convolved_fluence_map_dimension[0] - 4) / 2 * this->pixel_size, \
        (float)(this->convolved_fluence_map_dimension[1] - 4) / 2 * this->pixel_size, -this->SAD};
    float scaling_factor = (this->sampling_range[1] - this->sampling_range[0]) \
        / (this->sampling_points - 1) / this->SAD;
    vector_z_prime[0] = vector_z_prime[0] * scaling_factor;
    vector_z_prime[1] = vector_z_prime[1] * scaling_factor;
    vector_z_prime[2] = vector_z_prime[2] * scaling_factor;

    float BEV_cell_diagnal[3]{vector_x[0] + vector_y[0] - vector_z_prime[0], \
        vector_x[1] + vector_y[1] - vector_z_prime[1], \
        vector_x[2] + vector_y[2] - vector_z_prime[2]};
    
    float diagnal = sqrt(BEV_cell_diagnal[0] * BEV_cell_diagnal[0] + \
        BEV_cell_diagnal[1] * BEV_cell_diagnal[1] + \
        BEV_cell_diagnal[2] * BEV_cell_diagnal[2]);
    
    if (diagnal >= Phtm.voxelSize)
    {
        cout << "The diagnal of the BEV cell is larger than the phantom dimension, which means "
            "that more than one phantom grid point may be within a BEV grid, incompatible with "
            "gradient calculation" << endl;
        exit;
    }
}

extern "C"
void BEVDoseForward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    cudaSurfaceObject_t dose_surface, \
    float phantom_size[3], float phantom_iso[3], \
    cudaTextureObject_t phantom_texture, \
    float* convolved_fluence_map, \
    FCBBkernel* FCBB_kernel, \
    cudaStream_t stream=0);

void beam::BEV_dose_forward(phantom& Phtm, FCBBkernel* FCBB_kernel, cudaStream_t stream)
{
    float phantom_size[3]{Phtm.dimension[0] * Phtm.voxelSize, Phtm.dimension[1] * \
        Phtm.voxelSize, Phtm.pitch * Phtm.voxelSize};
    float phantom_iso[3]{Phtm.isocenter[0], Phtm.isocenter[1], Phtm.isocenter[2]};
    BEVDoseForward(this->zenith, this->azimuth, this->SAD, this->pixel_size, \
        this->sampling_range[0], this->sampling_range[1], this->sampling_points, \
        this->FCBB_BEV_dose_surface, \
        phantom_size, phantom_iso, \
        Phtm.tex, \
        this->d_convolved_fluence_map, \
        FCBB_kernel, \
        stream);
}

extern "C"
void writeSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);
extern "C"
void readTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* data);

void E2E::test_volume_rendering()
{
    // array<int, 3> volumeSize({128, 128, 128});
    array<int, 3> volumeSize({64, 64, 128});
    uint volume = volumeSize[0] * volumeSize[2] * volumeSize[2];
    float* h_volume = (float*)malloc(volume*sizeof(float));
    string volumeIn{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/volumeIn.dat"};
    ifstream inFile(volumeIn);
    inFile.read((char*)h_volume, volume*sizeof(float));
    inFile.close();

    float* d_volume;
    checkCudaErrors(cudaMalloc(&d_volume, volume*sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_volume, h_volume, volume*sizeof(float), cudaMemcpyHostToDevice));

    // volume_ initialization
    cudaExtent volumeSize_ = make_cudaExtent(volumeSize[2], volumeSize[1], volumeSize[0]);
    cudaArray* content;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    checkCudaErrors(cudaMalloc3DArray(&content, &channelDesc, volumeSize_, cudaArraySurfaceLoadStore));

    cudaSurfaceObject_t volumeSurf;
    cudaResourceDesc surfRes;
    memset(&surfRes, 0, sizeof(cudaResourceDesc));
    surfRes.resType = cudaResourceTypeArray;
    surfRes.res.array.array = content;
    checkCudaErrors(cudaCreateSurfaceObject(&volumeSurf, &surfRes));

    cudaTextureObject_t volumeTex;
    cudaResourceDesc texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array = content;

    cudaTextureDesc texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));
    texDescr.normalizedCoords = true;
    texDescr.filterMode = cudaFilterModeLinear;
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.addressMode[1] = cudaAddressModeBorder;
    texDescr.addressMode[2] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;

    checkCudaErrors(cudaCreateTextureObject(&volumeTex, &texRes, &texDescr, NULL));

    dim3 blockSize(1, 16, 16);
    dim3 gridSize(volumeSize[0] / blockSize.x, volumeSize[1] / blockSize.y, volumeSize[2] / blockSize.z);
    writeSurface(gridSize, blockSize, volumeSurf, d_volume);

    array<int, 3> outputVolume({96, 96, 96});
    volume = outputVolume[0] * outputVolume[1] * outputVolume[2];
    gridSize = dim3(outputVolume[0] / blockSize.x, outputVolume[1] / blockSize.y, \
        outputVolume[2] / blockSize.z);
    float* h_output = (float*)malloc(volume*sizeof(float));
    float* d_output;
    checkCudaErrors(cudaMalloc(&d_output, volume*sizeof(float)));
    readTexture(gridSize, blockSize, volumeTex, d_output);
    checkCudaErrors(cudaMemcpy(h_output, d_output, volume*sizeof(float), cudaMemcpyDeviceToHost));

    string volumeOut{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/volumeOut.dat"};
    ofstream outFile(volumeOut);
    outFile.write((char*)h_output, volume*sizeof(float));
    outFile.close();
}


extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::test_BEV_dose_forward()
{
    // phantom initialization
    phantom Phtm;
    phantom_init_default(Phtm);
    Phtm.to_device();
    Phtm.textureInit();

    // beam initialization
    beam Beam;
    Beam.zenith = PI / 2;
    Beam.azimuth = 0;
    Beam.SAD = get_args<float>("SAD") / 10;
    Beam.pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
    Beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
    Beam.convolved_fluence_map_dimension = array<int, 2>({FM_dimension + 2 * FM_convolution_radius, \
        FM_dimension + 2 * FM_convolution_radius});
    Beam.extended_fluence_map_dimension = array<int, 2>({FM_dimension + 4 * FM_convolution_radius, \
        FM_dimension + 4 * FM_convolution_radius});
    Beam.FCBBinit(Phtm);

    uint volume = Beam.convolved_fluence_map_dimension[0] * Beam.convolved_fluence_map_dimension[1];
    checkCudaErrors(cudaMalloc(&(Beam.d_convolved_fluence_map), volume * sizeof(float)));
    float* h_convolved_fluence_map = (float*)malloc(volume * sizeof(float));
    for (uint i=0; i< volume; i++)
        h_convolved_fluence_map[i] = 1;
    checkCudaErrors(cudaMemcpy(Beam.d_convolved_fluence_map, h_convolved_fluence_map, \
        volume*sizeof(float), cudaMemcpyHostToDevice));
    
    (*FCBB6MeV).texInit();
    
    // // for debug purposes
    // uint volume_debug = volume * Beam.sampling_points;
    // checkCudaErrors(cudaMalloc(&HU_debug, volume_debug*sizeof(float)));
    // checkCudaErrors(cudaMalloc(&dose_debug, volume_debug*sizeof(float)));

    Beam.BEV_dose_forward(Phtm);

    // // for debug purposes
    // float* h_HU_debug = (float*)malloc(volume_debug*sizeof(float));
    // float* h_dose_debug = (float*)malloc(volume_debug*sizeof(float));
    // checkCudaErrors(cudaMemcpy(h_HU_debug, HU_debug, volume_debug*sizeof(float), cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(h_dose_debug, dose_debug, volume_debug*sizeof(float), cudaMemcpyDeviceToHost));
    // string h_HU_debug_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/h_HU_debug.dat"};
    // string h_dose_debug_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/h_dose_debug.dat"};
    // ofstream outPath(h_HU_debug_path);
    // outPath.write((char*)h_HU_debug, volume_debug*sizeof(float));
    // outPath.close();
    // outPath.open(h_dose_debug_path);
    // outPath.write((char*)h_dose_debug, volume_debug*sizeof(float));
    // outPath.close();


    // logical order : (z, x, y)
    array<int, 3> content_dim({Beam.sampling_points, Beam.convolved_fluence_map_dimension[0], \
        Beam.convolved_fluence_map_dimension[1]});
    volume = content_dim[0] * content_dim[1] * content_dim[2];
    float* h_content = (float*)malloc(volume*sizeof(float));
    float* d_content = nullptr;
    checkCudaErrors(cudaMalloc(&d_content, volume*sizeof(float)));

    uint blockS = 16;
    dim3 blockSize(blockS, blockS, 1);
    dim3 gridSize(Beam.convolved_fluence_map_dimension[0] / blockS, \
        Beam.convolved_fluence_map_dimension[1] / blockS, Beam.sampling_points);
    readSurface(gridSize, blockSize, Beam.FCBB_BEV_dose_surface, d_content);
    checkCudaErrors(cudaMemcpy(h_content, d_content, volume*sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/BEVforward.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_content, volume*sizeof(float));
    outFile.close();
}