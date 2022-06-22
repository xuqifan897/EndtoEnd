#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"

using namespace E2E;
using namespace std;

extern "C" void
PVCSDoseForward(float voxel_size, uint phantom_dim[3], uint phantom_pitch, \
    float phantom_iso[3], \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_size_x, float fluence_map_size_y, 
    float* FCBB_PVCS_dose, \
    cudaTextureObject_t BEV_dose_texture,
    cudaStream_t stream=0);

void beam::PVCS_dose_forward(phantom& Phtm)
{
    uint phantom_dim[3]{(Phtm.dimension)[0], (Phtm.dimension)[1], (Phtm.dimension)[2]};
    float phantom_iso[3]{Phtm.isocenter[0], Phtm.isocenter[1], Phtm.isocenter[2]};
    float fluence_map_size[2]{this->convolved_fluence_map_dimension[0] * this->pixel_size, \
        this->convolved_fluence_map_dimension[1] * this->pixel_size};
    PVCSDoseForward(Phtm.voxelSize, phantom_dim, Phtm.pitch, \
        phantom_iso, \
        this->zenith, this->azimuth, this->SAD, \
        this->sampling_range[0], this->sampling_range[1], this->sampling_points, \
        fluence_map_size[0], fluence_map_size[1], \
        this->d_FCBB_PVCS_dose, \
        this->FCBB_BEV_dose_texture);
}

extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::test_PVCS_dose_forward()
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
    Beam.BEV_dose_forward(Phtm);


    // // for debug purposes
    // array<int, 3> content_dim({Beam.sampling_points, Beam.convolved_fluence_map_dimension[0], \
    //     Beam.convolved_fluence_map_dimension[1]});
    // volume = content_dim[0] * content_dim[1] * content_dim[2];
    // float* h_content = (float*)malloc(volume*sizeof(float));
    // float* d_content = nullptr;
    // checkCudaErrors(cudaMalloc(&d_content, volume*sizeof(float)));

    // uint blockS = 16;
    // dim3 blockSize(blockS, blockS, 1);
    // dim3 gridSize(Beam.convolved_fluence_map_dimension[0] / blockS, \
    //     Beam.convolved_fluence_map_dimension[1] / blockS, Beam.sampling_points);
    // readSurface(gridSize, blockSize, Beam.FCBB_BEV_dose_surface, d_content);
    // checkCudaErrors(cudaMemcpy(h_content, d_content, volume*sizeof(float), cudaMemcpyDeviceToHost));
    // string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/PVCS_BEV_debug.dat"};
    // ofstream outFile(outputPath);
    // outFile.write((char*)h_content, volume*sizeof(float));
    // outFile.close();

    Beam.PVCS_dose_forward(Phtm);

    volume = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_FCBB_PVCS_dose = (float*)malloc(volume * sizeof(float));
    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose, Beam.d_FCBB_PVCS_dose, \
        volume*sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/PVCSforward.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_FCBB_PVCS_dose, volume*sizeof(float));
    outFile.close();
}


extern "C"
void testWritePVCSSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface);

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output);


void E2E::test_PVCS_surface()
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

    dim3 blockSize(16, 16, 1);
    dim3 gridSize;
    gridSize.x = ceil((float)(Phtm.dimension[0]) / blockSize.x);
    gridSize.y = ceil((float)(Phtm.dimension[1]) / blockSize.y);
    gridSize.z = ceil((float)(Phtm.pitch) / blockSize.z);
    testWritePVCSSurface(gridSize, blockSize, Beam.FCBB_PVCS_dose_grad_surface);

    uint volume = blockSize.x * blockSize.y * blockSize.z * gridSize.x * gridSize.y * gridSize.z;
    float* h_output = (float*)malloc(volume*sizeof(float));
    float* d_output = nullptr;
    checkCudaErrors(cudaMalloc(&d_output, volume*sizeof(float)));
    testReadPVCSTexture(gridSize, blockSize, Beam.FCBB_PVCS_dose_grad_texture, d_output);
    checkCudaErrors(cudaMemcpy(h_output, d_output, volume*sizeof(float), cudaMemcpyDeviceToHost));

    string test_PVCS_output{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/test_PVCS_surface.dat"};
    ofstream outFile(test_PVCS_output);
    outFile.write((char*)h_output, volume*sizeof(float));
    outFile.close();
}