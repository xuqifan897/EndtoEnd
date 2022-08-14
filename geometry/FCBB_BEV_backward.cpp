#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "geom.h"
#include "args.h"
#include "optim.h"

using namespace E2E;
using namespace std;

extern "C"
void BEVDoseBackward(float zenith, float azimuth, float SAD, float pixel_size, \
    float sampling_range_start, float sampling_range_end, uint sampling_points, \
    float phantom_size[3], float phantom_iso[3], \
    float* d_convolved_fluence_grad, \
    cudaTextureObject_t phantom_texture, \
    float* d_FCBB_BEV_dose_grad, \
    FCBBkernel* FCBB_kernel, \
    cudaStream_t stream);

void beam::BEV_dose_backward(phantom& Phtm, FCBBkernel* kernel, cudaStream_t stream)
{
    float phantom_size[3]{Phtm.dimension[0]*Phtm.voxelSize, Phtm.dimension[1]*\
        Phtm.voxelSize, Phtm.pitch*Phtm.voxelSize};
    float phantom_iso[3]{this->isocenter[0], this->isocenter[1], this->isocenter[2]};
    if (this->d_convolved_fluence_map_grad == nullptr)
    {
        cout << "d_convolved_fluence_map_grad is not initialized. E2E::beams_init() not called." << endl;
        exit;
    }
    BEVDoseBackward(this->zenith, this->azimuth, this->SAD, this->pixel_size, \
        this->sampling_range[0], this->sampling_range[1], this->sampling_points, \
        phantom_size, phantom_iso, \
        this->d_convolved_fluence_map_grad, \
        Phtm.tex, \
        this->d_FCBB_BEV_dose_grad, \
        kernel, \
        stream);
}

extern "C"
void readSurface(dim3 gridSize, dim3 blockSize, cudaSurfaceObject_t surface, float* data);

void E2E::test_FCBB_BEV_backward(std::vector<beam>& beams, phantom& Phtm)
{
    if (!beam::FCBB_PVCS_dose_grad_init)
        beam::FCBBStaticInit(Phtm);
    
    beam& this_beam = beams[0];
    this_beam.FCBBinit(Phtm);

    uint part = 2;

    if (part == 1)
    {
        // first, we conduct random number experiment
        uint BEV_dose_grad_size = this_beam.convolved_fluence_map_dimension[0] * \
            this_beam.convolved_fluence_map_dimension[1] * this_beam.sampling_points;
        float* h_FCBB_BEV_dose_grad = (float*)malloc(BEV_dose_grad_size*sizeof(float));
        srand(1008611);
        uint norm = 1024;
        for (uint i=0; i<BEV_dose_grad_size; i++)
            h_FCBB_BEV_dose_grad[i] = (float)(rand()%norm) / (norm-1);
        checkCudaErrors(cudaMemcpy(this_beam.d_FCBB_BEV_dose_grad, h_FCBB_BEV_dose_grad, \
            BEV_dose_grad_size*sizeof(float), cudaMemcpyHostToDevice));
        
        this_beam.BEV_dose_backward(Phtm);

        // read result
        uint convolved_fluence_map_size = this_beam.convolved_fluence_map_dimension[0] * \
            this_beam.convolved_fluence_map_dimension[1];
        float* h_convolved_fluence_map_grad = (float*)malloc(convolved_fluence_map_size*sizeof(float));
        checkCudaErrors(cudaMemcpy(h_convolved_fluence_map_grad, this_beam.d_convolved_fluence_map_grad, \
            convolved_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
        string outfile{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/convolved_grad_rand.dat"};
        ofstream outFile(outfile);
        if (! outFile.is_open())
        {
            cout << "Could not open this file: " << outfile << endl;
            exit;
        }
        outFile.write((char*)h_convolved_fluence_map_grad, convolved_fluence_map_size*sizeof(float));
        outFile.close();
    }
    else if (part == 2)
    {
        // in the second part, we verify the correctness of BEV backward
        // in the first step, we initialize the convolved_fluence_map, then do forward
        uint convolved_fluence_map_size = this_beam.convolved_fluence_map_dimension[0] * \
            this_beam.convolved_fluence_map_dimension[1];
        float* h_convolved_fluence_map = (float*)malloc(convolved_fluence_map_size*sizeof(float));
        for (uint i=0; i< convolved_fluence_map_size; i++)
            h_convolved_fluence_map[i] = 1;
        checkCudaErrors(cudaMemcpy(this_beam.d_convolved_fluence_map, h_convolved_fluence_map, \
            convolved_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
        
        this_beam.BEV_dose_forward(Phtm);

        uint BEV_dose_size = convolved_fluence_map_size * this_beam.sampling_points;
        float* d_BEV_dose;
        float* h_BEV_dose = (float*)malloc(BEV_dose_size*sizeof(float));
        checkCudaErrors(cudaMalloc((void**)&d_BEV_dose, BEV_dose_size*sizeof(float)));

        dim3 blockSize(8, 8, 8);
        dim3 gridSize(this_beam.convolved_fluence_map_dimension[0] / blockSize.x, \
            this_beam.convolved_fluence_map_dimension[1] / blockSize.y, this_beam.sampling_points / blockSize.z);
        readSurface(gridSize, blockSize, this_beam.FCBB_BEV_dose_surface, d_BEV_dose);
        checkCudaErrors(cudaMemcpy(h_BEV_dose, d_BEV_dose, BEV_dose_size*sizeof(float), cudaMemcpyDeviceToHost));

        float* h_BEV_dose_grad = (float*)malloc(BEV_dose_size*sizeof(float));
        float* h_convolved_fluence_map_grad = (float*)malloc(convolved_fluence_map_size * sizeof(float));
        srand(1008611);
        uint num_samples = 100; // examine num_samples points
        for (uint i=0; i<num_samples; i++)
        {
            uint idx = rand() % BEV_dose_size;
            float forward_value = h_BEV_dose[idx];

            // clear h_BEV_dose_grad
            for (uint j=0; j<BEV_dose_size; j++)
                h_BEV_dose_grad[j] = 0;
            h_BEV_dose_grad[idx] = 1;
            checkCudaErrors(cudaMemcpy(this_beam.d_FCBB_BEV_dose_grad, h_BEV_dose_grad, \
                BEV_dose_size*sizeof(float), cudaMemcpyHostToDevice));

            this_beam.BEV_dose_backward(Phtm);

            checkCudaErrors(cudaMemcpy(h_convolved_fluence_map_grad, this_beam.d_convolved_fluence_map_grad, \
                convolved_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
            uint convolved_fluence_map_idx = idx % convolved_fluence_map_size;
            float backward_value = h_convolved_fluence_map_grad[convolved_fluence_map_idx];

            cout << "forward: " << forward_value << " " << "backward: " << backward_value << endl;
        }

        // string BEV_dose_path{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary/BEV_unit_forward.dat"};
        // ofstream outFile(BEV_dose_path);
        // if (! outFile.is_open())
        // {
        //     cout << "Could not open this file: " << BEV_dose_path << endl;
        //     exit;
        // }
        // outFile.write((char*)h_BEV_dose, BEV_dose_size*sizeof(float));
        // outFile.close();
    }
}