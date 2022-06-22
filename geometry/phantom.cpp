#include <array>
#include <fstream>
#include <iostream>
#include "geom.h"
#include "args.h"
#include <string>
#include <vector>
#include <math.h>

#include <cuda_runtime.h>
#include <helper_cuda.h>

using namespace E2E;
using namespace std;

phantom::phantom()
{   
    // this->dimension_org = array<int, 3>({0, 0, 0});
    this->dimension = array<int, 3>({0, 0, 0});
    this->isocenter = array<float, 3>({0, 0, 0});
    this->voxelSize = 0;

    this->pitch_module = E2E::pitch_module;
    this->pitch = 0;
    this->pitchPadding = false;

    this->h_HU = nullptr;
    this->h_PTVweight = nullptr;
    this->h_PTVtarget = nullptr;
    this->h_OARweight = nullptr;
    this->h_OARtarget = nullptr;
    this->h_Dose = nullptr;

    this->d_HU = 0;
    this->d_PTVweight = nullptr;
    this->d_PTVtarget = nullptr;
    this->d_OARweight = nullptr;
    this->d_OARtarget = nullptr;
    this->d_Dose = nullptr;

    this->texInit = false;
}

phantom::~phantom()
{
    if (this->h_HU != nullptr)
        free(this->h_HU);
    if (this->h_PTVweight != nullptr)
        free(this->h_PTVweight);
    if (this->h_PTVtarget != nullptr)
        free(this->h_PTVtarget);
    if (this->h_OARweight != nullptr)
        free(this->h_OARweight);
    if (this->h_OARtarget != nullptr)
        free(this->h_OARtarget);
    if (this->h_Dose != nullptr)
        free(this->h_Dose);
    
    textureDecon();
    
    if (this->d_HU != 0)
        checkCudaErrors(cudaFreeArray(this->d_HU));
    if (this->d_PTVweight != nullptr)
        checkCudaErrors(cudaFree(this->d_PTVweight));
    if (this->d_PTVtarget != nullptr)
        checkCudaErrors(cudaFree(this->d_PTVtarget));
    if (this->d_OARweight != nullptr)
        checkCudaErrors(cudaFree(this->d_OARweight));
    if (this->d_OARtarget != nullptr)
        checkCudaErrors(cudaFree(this->d_OARtarget));
    if (this->d_Dose != nullptr)
        checkCudaErrors(cudaFree(this->d_Dose));
}

void phantom::pitchPad()
{
    this->pitchPadding = true;
    this->pitch = ceil((float)(this->dimension[2]) / this->pitch_module) * this->pitch_module;
    size_t new_size = this->dimension[0] * this->dimension[1] * this->pitch;
    vector<float**> pointers{&(this->h_HU), &(this->h_PTVweight), &(this->h_PTVtarget), &(this->h_OARweight), &(this->h_OARtarget)};
    for (int i=0; i<pointers.size(); i++)
    {
        if (*(pointers[i]) == nullptr)
            continue;
        
        float* ptr_temp = (float*)malloc(new_size * sizeof(float));
        for (int j=0; j<this->dimension[0]; j++)
        {   
            size_t jXdimension1 = j * this->dimension[1];
            for (int k=0; k<this->dimension[1]; k++)
            {
                size_t Xdimension2 = (jXdimension1 + k) * this->dimension[2];
                size_t Ydimension2 = (jXdimension1 + k) * this->pitch;
                for (int l=0; l<this->dimension[2]; l++)
                    ptr_temp[Ydimension2 + l] = (*(pointers[i]))[Xdimension2 + l];
            }
        }

        free(*(pointers[i]));
        *(pointers[i]) = ptr_temp;
    }
}

// void phantom::pitchPad()
// {
//     this->pitchPadding = true;
//     array<int, 3> new_dimension;
//     new_dimension[0] = ceil((float)(this->dimension[0]) / this->pitch_module) * this->pitch_module;
//     new_dimension[1] = ceil((float)(this->dimension[1]) / this->pitch_module) * this->pitch_module;
//     new_dimension[2] = ceil((float)(this->dimension[2]) / this->pitch_module) * this->pitch_module;
//     this->pitch = new_dimension[2];

//     size_t new_size = new_dimension[0] * new_dimension[1] * new_dimension[2];
//     vector<float**> pointers{&(this->h_HU), &(this->h_PTVweight), &(this->h_PTVtarget), &(this->h_OARweight), &(this->h_OARtarget)};
//     for (int i=0; i<pointers.size(); i++)
//     {
//         if (*(pointers[i]) == nullptr)
//             continue;
        
//         float* ptr_temp = (float*)malloc(new_size * sizeof(float));
//         for (int j=0; j<new_size; j++)
//             ptr_temp[j] = 0;
        
//         for (int j=0; j<this->dimension[0]; j++)
//         {
//             size_t Xterm = j * this->dimension[1];
//             size_t new_Xterm = j * new_dimension[1];
//             for (int k=0; k<this->dimension[1]; k++)
//             {
//                 size_t Yterm = (Xterm + k) * this->dimension[2];
//                 size_t new_Yterm = (new_Xterm + k) * new_dimension[2];
//                 for (int l=0; l<this->dimension[2]; l++)
//                 {
//                     size_t idx = Yterm + l;
//                     size_t new_idx = new_Yterm + l;
//                     ptr_temp[new_idx] = (*(pointers[i]))[idx];
//                 }
//             }
//         }
//         free(*(pointers[i]));
//         *(pointers[i]) = ptr_temp;
//     }
//     this->dimension_org = this->dimension;
//     this->dimension = new_dimension;
// }

void phantom::Dose_init()
{
    if (! this->pitchPadding)
    {
        cout << "pitchPad() member function has not been called" << endl;
        exit;
    }
    size_t size = this->dimension[0] * this->dimension[1] * this->pitch;
    this->h_Dose = (float*)malloc(size*sizeof(float));
    for (int i=0; i<size; i++)
        this->h_Dose[i] = 0;
}

void phantom::to_device()
{
    if (! this->pitchPadding)
    {
        cout << "pitchPad() member function has not been called" << endl;
        exit;
    }
    vector<float**> h_pointers{&(this->h_PTVweight), &(this->h_PTVtarget),\
        &(this->h_OARweight), &(this->h_OARtarget), &(this->h_Dose)};
    vector<float**> d_pointers{&(this->d_PTVweight), &(this->d_PTVtarget),\
        &(this->d_OARweight), &(this->d_OARtarget), &(this->d_Dose)};
    size_t size = this->dimension[0] * this->dimension[1] * this->pitch * sizeof(float);
    for(int i=0; i<5; i++)
    {
        if (*(d_pointers[i]) == nullptr)
            checkCudaErrors(cudaMalloc(d_pointers[i], size));
        checkCudaErrors(cudaMemcpy(*(d_pointers[i]), *(h_pointers[i]),\
            size, cudaMemcpyHostToDevice));
    }

    const cudaExtent volumeSize_ = make_cudaExtent(\
        this->pitch, this->dimension[1], this->dimension[0]);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    checkCudaErrors(cudaMalloc3DArray(&(this->d_HU), &channelDesc, volumeSize_));

    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr((void*)(this->h_HU),\
        volumeSize_.width*sizeof(float), volumeSize_.width, (this->dimension)[1]);
    copyParams.dstArray = this->d_HU;
    copyParams.extent = volumeSize_;
    copyParams.kind = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));
}

void phantom::DoseToHost()
{
    size_t size = this->dimension[0] * this->dimension[1] * this->pitch * sizeof(float);
    if (this->h_Dose == nullptr)
        this->h_Dose = (float*)malloc(size);
    checkCudaErrors(cudaMemcpy(this->h_Dose, this->d_Dose, size, cudaMemcpyDeviceToHost));
}

void phantom::textureInit()
{
    this->texInit = true;
    cudaResourceDesc texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array = this->d_HU;

    cudaTextureDesc texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));

    texDescr.normalizedCoords = true; // access with normalized texture coordinates
    texDescr.filterMode = cudaFilterModeLinear; // linear interpolation
    // wrap texture coordinates
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.addressMode[1] = cudaAddressModeBorder;
    texDescr.addressMode[2] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;

    checkCudaErrors(cudaCreateTextureObject(&(this->tex), &texRes, &texDescr, NULL));
}

void phantom::textureDecon()
{
    if (this->tex)
        checkCudaErrors(cudaDestroyTextureObject(this->tex));
}

int E2E::phantom_init_default(phantom& Phtm)
{
    vector<int> dim = get_args<vector<int>>("phantom-dimension");
    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    float voxelsize = get_args<float>("voxel-size");
    string phantom_path = get_args<string>("phantom-path");
    string PTV_weight_path = get_args<string>("PTV-weight-path");
    string PTV_target_path = get_args<string>("PTV-target-path");
    string OAR_weight_path = get_args<string>("OAR-weight-path");
    string OAR_target_path = get_args<string>("OAR-target-path");

    Phtm.dimension = array<int, 3>({dim[0], dim[1], dim[2]});
    Phtm.isocenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});
    Phtm.voxelSize = voxelsize/10;
    int size = dim[0] * dim[1] * dim[2];

    vector<float**> pointers{&(Phtm.h_HU), &(Phtm.h_PTVweight), &(Phtm.h_PTVtarget), \
        &(Phtm.h_OARweight), &(Phtm.h_OARtarget)};
    vector<string> files{phantom_path, PTV_weight_path, PTV_target_path, \
        OAR_weight_path, OAR_target_path};
    for (int i=0; i<5; i++)
    {
        *(pointers[i]) = (float*)malloc(size*sizeof(float));
        ifstream input_file(files[i]);
        if (! input_file.is_open())
        {
            cerr << "Could not open this file: " << files[i] << endl;
            return EXIT_FAILURE;
        }
        input_file.read((char*)(*(pointers[i])), size*sizeof(float));
    }

    // normalize
    for (int i=0; i<size; i++)
        Phtm.h_HU[i] /= 1000;
    
    Phtm.pitchPad();
    Phtm.Dose_init();
    return 0;
}

extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, \
    float* output, int y_points, int z_points, float y_scale, float z_scale, \
    float x, cudaTextureObject_t& texObj);

void E2E::runTest(phantom& Phtm)
{
    vector<int> dim = get_args<vector<int>>("phantom-dimension"); // (x, y, z)
    int y_points = ceil((float)dim[1] / Phtm.pitch_module) * Phtm.pitch_module;
    int z_points = ceil((float)dim[2] / Phtm.pitch_module) * Phtm.pitch_module;
    float x = 0.7;
    
    float* h_slice = (float*)malloc(y_points*z_points*sizeof(float));
    float* d_slice = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_slice, y_points*z_points*sizeof(float)));

    const dim3 blockSize(Phtm.pitch_module, Phtm.pitch_module);
    const dim3 gridSize(y_points/Phtm.pitch_module, z_points/Phtm.pitch_module);
    render_kernel(gridSize, blockSize, d_slice, y_points, z_points, \
        1, (float)Phtm.dimension[2] / Phtm.pitch, x, Phtm.tex);

    checkCudaErrors(cudaMemcpy(h_slice, d_slice, y_points*z_points*sizeof(float),\
        cudaMemcpyDeviceToHost));
    string filename {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/slice01.dat"};
    ofstream outFile(filename);
    if (! outFile.is_open())
    {
        cout << "Could not open this file: " << filename << endl;
        exit;
    }
    outFile.write((char*)h_slice, y_points*z_points*sizeof(float));
    outFile.close();
}