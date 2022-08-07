#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <fstream>
#include <vector>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

extern "C"
void FCBBPVCSDoseGrad(cudaSurfaceObject_t FCBB_PVCS_dose_grad_surface, float* elementWiseLoss, \
    float* d_FCBB_PVCS_dose, float* PTV_weight, float* PTV_target, \
    float* OAR_weight, float* OAR_target, uint dimension[3], cudaStream_t stream);

void beam::calc_FCBB_PVCS_dose_grad(phantom& Phtm, float** d_elementWiseLoss, \
        float* d_PVCS_total_dose, cudaStream_t stream)
{
    /* the loss function is defined as $$PTV_weight .* ||PVCS_dose - PTV_target||_2^2 +
     OAR_weight .* ||(PVCS_dose - OAR_target)_+||_2^2$$. So the loss function should be 
     2 * PTV_weight .* (PVCS_dose - PTV_target) + OAR_weight .* (PVCS_dose - OAR_target)_+ */
    
    if (! Phtm.pitchPadding)
    {
        cout << "It is required that pitchPad() must be called" << endl;
        exit;
    }

    if (! FCBB_PVCS_dose_grad_init)
    {
        cout << "FCBBStaticInit static member function is not called!" << endl;
        exit;
    }

    if (*d_elementWiseLoss == nullptr)
    {
        cout << "d_elementWiseLoss has not been initialized" << endl;
        exit;
    }
    uint dimension[3]{Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch};
    FCBBPVCSDoseGrad(FCBB_PVCS_dose_grad_surface, *d_elementWiseLoss, \
        d_PVCS_total_dose, Phtm.d_PTVweight, Phtm.d_PTVtarget, \
        Phtm.d_OARweight, Phtm.d_OARtarget, dimension, stream);
}

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output);

void E2E::test_calc_FCBB_PVCS_dose_grad(vector<beam>& beams, phantom& Phtm)
{
    if (! beam::FCBB_PVCS_dose_grad_init)
        beam::FCBBStaticInit(Phtm);

    beams[0].FCBBinit(Phtm);

    string inputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/doseRand.dat"};
    uint size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* rand = (float*)malloc(size*5*sizeof(float));
    ifstream inFile(inputPath);
    inFile.read((char*)rand, size*5*sizeof(float));
    inFile.close();

    float* h_FCBB_PVCS_dose = rand;
    float* h_PTVweight = rand + size;
    float* h_PTVtarget = h_PTVweight + size;
    float* h_OARweight = h_PTVtarget + size;
    float* h_OARtarget = h_OARweight + size;

    checkCudaErrors(cudaMemcpy(beams[0].d_FCBB_PVCS_dose, h_FCBB_PVCS_dose, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVweight, h_PTVweight, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_PTVtarget, h_PTVtarget, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARweight, h_OARweight, size*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(Phtm.d_OARtarget, h_OARtarget, size*sizeof(float), cudaMemcpyHostToDevice));

    float* d_elementWiseLoss = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_elementWiseLoss, size*sizeof(float)));
    float* h_elementWiseLoss = (float*)malloc(size*sizeof(float));

    // // for debug purposes
    // checkCudaErrors(cudaMalloc((void**)&dose_debug, size*sizeof(float)));

    beam::calc_FCBB_PVCS_dose_grad(Phtm, &d_elementWiseLoss, beams[0].d_FCBB_PVCS_dose);
    // float reduce_value = reduction(d_elementWiseLoss, size);
    // cout << reduce_value << endl;

    // // for debug purposes
    // float* h_dose_debug = (float*)malloc(size*sizeof(float));
    // checkCudaErrors(cudaMemcpy(h_dose_debug, dose_debug, size*sizeof(float), cudaMemcpyDeviceToHost));
    // string outputPath_{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_dose_grad_.dat"};
    // ofstream outFile_(outputPath_);
    // outFile_.write((char*)h_dose_debug, size*sizeof(float));
    // outFile_.close();

    float* h_FCBB_PVCS_dose_grad = (float*)malloc(size*sizeof(float));
    float* d_FCBB_PVCS_dose_grad = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_PVCS_dose_grad, size*sizeof(float)));
    
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(Phtm.dimension[0] / blockSize.x, Phtm.dimension[1] / blockSize.y, \
        Phtm.pitch / blockSize.z);
    testReadPVCSTexture(gridSize, blockSize, beam::FCBB_PVCS_dose_grad_texture, d_FCBB_PVCS_dose_grad);

    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose_grad, d_FCBB_PVCS_dose_grad, size*sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_dose_grad.dat"};
    ofstream outFile(outputPath);
    outFile.write((char*)h_FCBB_PVCS_dose_grad, size*sizeof(float));
    outFile.close();

    checkCudaErrors(cudaMemcpy(h_elementWiseLoss, d_elementWiseLoss, size*sizeof(float), cudaMemcpyDeviceToHost));
    outputPath = "/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/FCBB_PVCS_element_wise_loss.dat";
    outFile.open(outputPath);
    outFile.write((char*)h_elementWiseLoss, size*sizeof(float));
    outFile.close();
}

extern "C"
void PVCSDoseBackward(float voxel_size, float phantom_iso[3], \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_pixel_size, uint fluence_map_dimension, \
    float* d_FCBB_BEV_dose_grad, cudaTextureObject_t PVCSDoseGradTexture, \
    cudaStream_t stream);

extern "C"
void PVCSDoseBackward_new(float voxel_size, float phantom_iso[3], \
    float zenith, float azimuth, float SAD, \
    float sampling_start, float sampling_end, uint sampling_points, \
    float fluence_map_pixel_size, uint fluence_map_dimension, \
    float* d_FCBB_BEV_dose_grad, cudaTextureObject_t PVCSDoseGradTexture, \
    cudaStream_t stream);

void beam::PVCS_dose_backward(phantom& Phtm, cudaStream_t stream)
{
    if (this->d_FCBB_BEV_dose_grad == nullptr)
    {
        cout << "d_FCBB_BEV_dose_grad is not initialized, beam::FCBBinit() not called" << endl;
        exit;
    }
    float phantom_iso[3]{this->isocenter[0], this->isocenter[1], this->isocenter[2]};
    PVCSDoseBackward_new(Phtm.voxelSize, phantom_iso, \
        this->zenith, this->azimuth, this->SAD, \
        this->sampling_range[0], this->sampling_range[1], this->sampling_points, \
        this->pixel_size, this->convolved_fluence_map_dimension[0], \
        this->d_FCBB_BEV_dose_grad, this->FCBB_PVCS_dose_grad_texture, \
        stream);
}

extern "C"
void writeFCBBPVCSDoseGradSurface(cudaSurfaceObject_t surface, float* input, \
    uint dim_x, uint dim_y, uint dim_z, cudaStream_t stream=0);

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output);

inline uint
get_idx(uint idx_x, uint idx_y, uint idx_z, uint dim_x, uint dim_y, uint dim_z)
{
    return (idx_x * dim_y + idx_y) * dim_z + idx_z;
}

inline void
get_coords(array<uint, 3>& coords, uint idx, uint dim_x, uint dim_y, uint dim_z)
{
    coords[0] = idx / (dim_y * dim_z);
    idx -= coords[0] * dim_y * dim_z;
    coords[1] = idx / dim_z;
    coords[2] = idx - coords[1] * dim_z;
}

void E2E::test_FCBB_PVCS_backward(std::vector<beam>& beams, phantom& Phtm)
{
    if (! beam::FCBB_PVCS_dose_grad_init)
        beam::FCBBStaticInit(Phtm);

    beams[0].FCBBinit(Phtm);

    uint part = 2;

    if (part == 1)
    {
        // first, we conduct trivial backward experiment. In other words, 
        // we set beam::FCBB_PVCS_dose_grad_surface to all 1
        uint input_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
        float* h_input = (float*)malloc(input_size * sizeof(float));
        for (uint i=0; i<input_size; i++)
            h_input[i] = 1;

        float* d_input = nullptr;
        checkCudaErrors(cudaMalloc((void**)&d_input, input_size*sizeof(float)));
        checkCudaErrors(cudaMemcpy(d_input, h_input, input_size*sizeof(float), cudaMemcpyHostToDevice));
        writeFCBBPVCSDoseGradSurface(beam::FCBB_PVCS_dose_grad_surface, d_input, \
            Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch, 0);
        
        beams[0].PVCS_dose_backward(Phtm);
        
        uint output_size = beams[0].convolved_fluence_map_dimension[0] * \
            beams[0].convolved_fluence_map_dimension[1] * beams[0].sampling_points;
        float* h_output = (float*)malloc(output_size*sizeof(float));
        checkCudaErrors(cudaMemcpy(h_output, beams[0].d_FCBB_BEV_dose_grad, \
            output_size*sizeof(float), cudaMemcpyDeviceToHost));
        
        string output_file{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/testFCBBPVCSBackward.dat"};
        ofstream outFile(output_file);
        if (! outFile.is_open())
        {
            cout << "Could not open this file: " << output_file << endl;
            exit;
        }
        outFile.write((char*)h_output, output_size*sizeof(float));
        outFile.close();

        // clean up
        free(h_input);
        free(h_output);
        checkCudaErrors(cudaFree(d_input));
    }
    else if(part == 2)
    {
        /* In this part, we randomly select some points in FCBB_BEV_dose_surface to 1, 
        then calculate its dose deposition to d_FCBB_PVCS_dose. In theory, the non-zero 
        values in d_FCBB_PVCS_dose should be the coefficients of trililnear interpolation. 
        So we iteratively set the corresponding elements in FCBB_PVCS_dose_grad_surface to 1, 
        then calculate the backward pass. The expected value of should be the coefficient.*/

        // firstly, we randomly selects points from PVCS space. i.e., x, y, z, range 1 .. 199
        beam& this_beam = beams[0];
        vector<array<uint, 3>> points{{6, 80, 120}, {32, 79, 110}, {50, 62, 95}, \
            {80, 140, 99}, {120, 128, 130}, {150, 97, 105}};
        uint BEV_size = this_beam.sampling_points * this_beam.convolved_fluence_map_dimension[0] * \
            this_beam.convolved_fluence_map_dimension[1];
        uint PVCS_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
        float* h_BEV_array = (float*)malloc(BEV_size * sizeof(float));
        float* h_PVCS_array = (float*)malloc(PVCS_size * sizeof(float));
        float* d_BEV_array = nullptr;
        checkCudaErrors(cudaMalloc((void**)&d_BEV_array, BEV_size*sizeof(float)));
        float* d_PVCS_array = nullptr;
        checkCudaErrors(cudaMalloc((void**)&d_PVCS_array, PVCS_size*sizeof(float)));

        for (uint i=0; i<points.size(); i++)
        {
            array<uint, 3>& point = points[i];
            uint temp_PVCS_idx = get_idx(point[0], point[1], point[2], \
                Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch);
            // clean up h_PVCS_array
            for (uint j=0; j< PVCS_size; j++)
                h_PVCS_array[j] = 0;
            h_PVCS_array[temp_PVCS_idx] = 1.;
            checkCudaErrors(cudaMemcpy(d_PVCS_array, h_PVCS_array, PVCS_size*sizeof(float), cudaMemcpyHostToDevice));
            writeFCBBPVCSDoseGradSurface(beam::FCBB_PVCS_dose_grad_surface, d_PVCS_array, \
                Phtm.dimension[0], Phtm.dimension[1], Phtm.pitch);
            
            this_beam.PVCS_dose_backward(Phtm);
            checkCudaErrors(cudaMemcpy(h_BEV_array, this_beam.d_FCBB_BEV_dose_grad, \
                BEV_size*sizeof(float), cudaMemcpyDeviceToHost));
            
            // theoretically, a PVCS point should have 8 associative BEV points (or less)
            vector<array<uint, 3>> non_zero_BEV_coords;
            vector<uint> non_zero_BEV_index;
            vector<float> non_zero_coefficients;
            for (uint j=0; j<BEV_size; j++)
            {
                if (h_BEV_array[j] != 0)
                {
                    non_zero_BEV_index.push_back(j);
                    non_zero_coefficients.push_back(h_BEV_array[j]);
                    non_zero_BEV_coords.push_back(array<uint, 3>());
                    get_coords(non_zero_BEV_coords.back(), j, this_beam.sampling_points, \
                        this_beam.convolved_fluence_map_dimension[0], this_beam.convolved_fluence_map_dimension[1]);
                    // cout << "(" << non_zero_BEV_coords.back()[0] << ", " << non_zero_BEV_coords.back()[1] << ", " \
                    //     << non_zero_BEV_coords.back()[2] << "): " << h_BEV_array[j] << endl;
                }
            }
            
            vector<float> test_coefficients;
            float std = 0;
            for (uint j=0; j<8; j++)
            {
                // clean up h_BEV_array
                for (uint k=0; k<BEV_size; k++)
                    h_BEV_array[k] = 0;
                h_BEV_array[non_zero_BEV_index[j]] = 1;
                checkCudaErrors(cudaMemcpy(d_BEV_array, h_BEV_array, BEV_size*sizeof(float), cudaMemcpyHostToDevice));
                writeFCBBPVCSDoseGradSurface(this_beam.FCBB_BEV_dose_surface, d_BEV_array, this_beam.sampling_points, \
                    this_beam.convolved_fluence_map_dimension[0], this_beam.convolved_fluence_map_dimension[1]);
                
                this_beam.PVCS_dose_forward(Phtm);

                checkCudaErrors(cudaMemcpy(h_PVCS_array, this_beam.d_FCBB_PVCS_dose, \
                    PVCS_size*sizeof(float), cudaMemcpyDeviceToHost));
                test_coefficients.push_back(h_PVCS_array[temp_PVCS_idx]);

                std += (non_zero_coefficients[j] - test_coefficients[j]) * (non_zero_coefficients[j] - test_coefficients[j]);
                cout << "(" << non_zero_BEV_coords[j][0] << ", " << non_zero_BEV_coords[j][1] << ", " \
                    << non_zero_BEV_coords[j][2] << "): " << non_zero_coefficients[j] << " " << \
                    test_coefficients[j] << endl;
            }
            std = sqrt(std / 8);
            cout << "standard deviation: " << std << endl << endl;
        }
    }
}

extern "C"
void textureMinusCoordinate(cudaTextureObject_t texture, float* d_output, \
    float coord0, float coord1, float coord2);
extern "C"
void writeSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::test_minus_coordinates_of_texture_memory_out_of_curiosity()
{
    array<int, 3> volumeSize({64, 64, 64});
    uint volume = volumeSize[0] * volumeSize[1] * volumeSize[2];
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

    float* d_output = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_output, sizeof(float)));
    float coord0 = - 1. / volumeSize[0] / 8;
    float coord1 = 0.5;
    float coord2 = 0.5;
    textureMinusCoordinate(volumeTex, d_output, coord0, coord1, coord2);

    float h_output;
    checkCudaErrors(cudaMemcpy(&h_output, d_output, sizeof(float), cudaMemcpyDeviceToHost));
    cout << h_output << endl;

    coord0 = 0.;
    textureMinusCoordinate(volumeTex, d_output, coord0, coord1, coord2);
    checkCudaErrors(cudaMemcpy(&h_output, d_output, sizeof(float), cudaMemcpyDeviceToHost));
    cout << h_output << endl;
}