#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

void fluence_map_init(vector<beam>& beams, phantom& Phtm)
{
    array<int, 2> extended_fluence_map_dimension{FM_dimension + 4 * FM_convolution_radius, \
        FM_dimension + 4 * FM_convolution_radius};
    uint extended_size = extended_fluence_map_dimension[0] * extended_fluence_map_dimension[1];
    // assign the initial value of the convolved fluence map to be 1/beams.size()
    float value = 1. / beams.size();
    float* h_extended_fluence_map = (float*)malloc(extended_size*sizeof(float));
    for (uint i=0; i<extended_size; i++)
        h_extended_fluence_map[i] = 0;
    for (uint i=2*FM_convolution_radius; i<extended_fluence_map_dimension[0] \
        -2*FM_convolution_radius; i++)
    {
        uint IDX = i * extended_fluence_map_dimension[1];
        for (uint j=2*FM_convolution_radius; j<extended_fluence_map_dimension[1] \
            -2*FM_convolution_radius; j++)
        {
            uint idx = IDX + j;
            h_extended_fluence_map[idx] = value;
        }
    }

    for (auto& beam : beams)
    {
        beam.FCBBinit(Phtm);
        checkCudaErrors(cudaMemcpy(beam.d_extended_fluence_map, h_extended_fluence_map, \
            extended_size*sizeof(float), cudaMemcpyHostToDevice));
    }
}

extern "C"
void testReadPVCSTexture(dim3 gridSize, dim3 blockSize, cudaTextureObject_t texture, float* output);

void optimize_stationary_view_results(vector<beam>& beams, phantom& Phtm, \
    float* d_PVCS_total_dose, float* d_element_wise_loss)
{
    // this function examines the forward results
    string results_dir{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary"};
    uint dose_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_PVCS_dose = (float*)malloc(dose_size * sizeof(float));
    for (int i=0; i<beams.size(); i++)
    {
        checkCudaErrors(cudaMemcpy(h_PVCS_dose, (beams[i].d_FCBB_PVCS_dose), \
            dose_size*sizeof(float), cudaMemcpyDeviceToHost));
        stringstream output_file_stream;
        output_file_stream << results_dir << "/dose" << i << ".dat";
        string output_file = output_file_stream.str();
        ofstream outFile(output_file);
        if (! outFile.is_open())
        {
            cout << "The file " << output_file << " could not open." << endl;
            exit;
        }
        outFile.write((char*)h_PVCS_dose, dose_size*sizeof(float));
        outFile.close();
    }

    checkCudaErrors(cudaMemcpy(h_PVCS_dose, d_PVCS_total_dose, \
        dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    stringstream output_file_stream;
    output_file_stream << results_dir << "/totalDose.dat";
    string output_file = output_file_stream.str();
    ofstream outFile(output_file);
    outFile.write((char*)h_PVCS_dose, dose_size*sizeof(float));
    outFile.close();

    checkCudaErrors(cudaMemcpy(h_PVCS_dose, d_element_wise_loss, \
        dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    output_file_stream.str(string());
    output_file_stream << results_dir << "/elementWiseLoss.dat";
    output_file = output_file_stream.str();
    outFile.open(output_file);
    outFile.write((char*)h_PVCS_dose, dose_size*sizeof(float));
    outFile.close();

    float* d_FCBB_PVCS_dose_grad = nullptr;
    float* h_FCBB_PVCS_dose_grad = (float*)malloc(dose_size*sizeof(float));
    checkCudaErrors(cudaMalloc((void**)&d_FCBB_PVCS_dose_grad, dose_size*sizeof(float)));
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(Phtm.dimension[0] / blockSize.x, Phtm.dimension[1] / blockSize.y, \
        Phtm.pitch / blockSize.z);
    testReadPVCSTexture(gridSize, blockSize, beam::FCBB_PVCS_dose_grad_texture, d_FCBB_PVCS_dose_grad);
    checkCudaErrors(cudaMemcpy(h_FCBB_PVCS_dose_grad, d_FCBB_PVCS_dose_grad, \
        dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    output_file_stream.str(string());
    output_file_stream << results_dir << "/FCBBPVCSDoseGrad.dat";
    output_file = output_file_stream.str();
    outFile.open(output_file);
    outFile.write((char*)h_FCBB_PVCS_dose_grad, dose_size*sizeof(float));
    outFile.close();
}

void optimize_stationary_view_backward_results(vector<beam>& beams)
{
    // this function examines the backward results
    uint convolved_fluence_map_size = beams[0].convolved_fluence_map_dimension[0] * \
        beams[0].convolved_fluence_map_dimension[1];
    uint fluence_map_size = beams[0].fluence_map_dimension[0] * beams[0].fluence_map_dimension[0];
    uint BEV_dose_size = convolved_fluence_map_size * beams[0].sampling_points;
    string results_dir{"/data/qifan/projects_qlyu/EndtoEnd3/data/optimize_stationary"};
    ofstream outFile;

    float* h_convolved_fluence_map_grad = (float*)malloc(convolved_fluence_map_size*sizeof(float));
    float* h_fluence_map_grad = (float*)malloc(fluence_map_size*sizeof(float));
    float* h_BEV_dose_grad = (float*)malloc(BEV_dose_size*sizeof(float));
    for (uint i=0; i<beams.size(); i++)
    {
        beam& this_beam = beams[i];
        checkCudaErrors(cudaMemcpy(h_convolved_fluence_map_grad, this_beam.d_convolved_fluence_map_grad, \
            convolved_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(h_fluence_map_grad, this_beam.d_fluence_grad, \
            fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(h_BEV_dose_grad, this_beam.d_FCBB_BEV_dose_grad, \
            BEV_dose_size*sizeof(float), cudaMemcpyDeviceToHost));
        stringstream output_file_stream;
        output_file_stream << results_dir << "/convolved_fluence_grad_" << i << ".dat";
        string output_file{output_file_stream.str()};
        outFile.open(output_file);
        if (! outFile.is_open())
        {
            cout << "Could not open this file: " << output_file << endl;
            exit;
        }
        outFile.write((char*)h_convolved_fluence_map_grad, convolved_fluence_map_size*sizeof(float));
        outFile.close();

        output_file_stream.str(string());
        output_file_stream << results_dir << "/fluence_grad_" << i << ".dat";
        output_file = output_file_stream.str();
        outFile.open(output_file);
        if (! outFile.is_open())
        {
            cout << "Could not open this file: " << output_file << endl;
            exit;
        }
        outFile.write((char*)h_fluence_map_grad, fluence_map_size*sizeof(float));
        outFile.close();

        output_file_stream.str(string());
        output_file_stream << results_dir << "/BEV_dose_grad_" << i << ".dat";
        output_file = output_file_stream.str();
        outFile.open(output_file);
        if (! outFile.is_open())
        {
            cout << "Could not open this file: " << output_file << endl;
            exit;
        }
        outFile.write((char*)h_BEV_dose_grad, BEV_dose_size*sizeof(float));
        outFile.close();
    }
}

void E2E::optimize_stationary(vector<beam>& beams, phantom& Phtm)
{
    // initialize the extended fluence map of beams
    // the valid fluence map zone is set to 1
    beam::FCBBStaticInit(Phtm);
    fluence_map_init(beams, Phtm);

    // intermediate data initialization
    // int iterations = get_args<float>("iterations");
    int iterations = 1;
    int iter = 0;
    uint phantom_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* d_PVCS_total_dose = nullptr;
    float* d_element_wise_loss = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_PVCS_total_dose, phantom_size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&d_element_wise_loss, phantom_size*sizeof(float)));
    float* d_out0 = nullptr; // intermediate result for reduction
    float* loss = nullptr; // the final loss array
    float* h_loss = (float*)malloc(iterations*sizeof(float));
    checkCudaErrors(cudaMalloc((void**)&d_out0, REDUCTION_BLOCK_SIZE*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&loss, iterations*sizeof(float)));

    // d_sources and h_sources are to the arrays to contain the d_FCBB_PVCS_dose
    float** h_sources = (float**)malloc(beams.size()*sizeof(float*));
    float** d_sources = nullptr;
    for (int i=0; i<beams.size(); i++)
    {
        if (beams[i].d_FCBB_PVCS_dose == nullptr)
        {
            cout << "The " << i << "th beam has not called FCBBinit()" << endl;
            exit;
        }
        h_sources[i] = beams[i].d_FCBB_PVCS_dose;
    }
    checkCudaErrors(cudaMalloc((void***)&d_sources, beams.size()*sizeof(float*)));
    checkCudaErrors(cudaMemcpy(d_sources, h_sources, beams.size()*sizeof(float*), \
        cudaMemcpyHostToDevice));

    // initialize streams and graphs
    vector<cudaStream_t> stream_beams(beams.size()+1); // the last stream for reduction
    vector<cudaEvent_t> event_beams(beams.size());
    vector<cudaEvent_t> new_event_beams(beams.size()+1);
    cudaGraph_t graph_beams;
    cudaGraphExec_t graphExec_beams;
    bool graph_beams_initialized = false;
    for (int i=0; i<stream_beams.size(); i++)
        checkCudaErrors(cudaStreamCreate(&(stream_beams[i])));
    for (int i=0; i<event_beams.size(); i++)
        checkCudaErrors(cudaEventCreate(&(event_beams[i])));
    for (int i=0; i<new_event_beams.size(); i++)
        checkCudaErrors(cudaEventCreate(&(new_event_beams[i])));
    
    // if graph_beams is not initialized, initialize it
    if (! graph_beams_initialized)
    {
        checkCudaErrors(cudaStreamBeginCapture(stream_beams[0], \
            cudaStreamCaptureModeGlobal));
        checkCudaErrors(cudaEventRecord(event_beams[0], stream_beams[0]));
        for (int i=1; i<beams.size(); i++)
            checkCudaErrors(cudaStreamWaitEvent(stream_beams[i], event_beams[0], 0));
        
        // concurrently calculate all beams
        for (int i=0; i<beams.size(); i++)
        {
            beams[i].convolve(FCBB6MeV, stream_beams[i]);
            beams[i].BEV_dose_forward(Phtm, FCBB6MeV, stream_beams[i]);
            beams[i].PVCS_dose_forward(Phtm, stream_beams[i]);
            if (i != 0)
            {
                checkCudaErrors(cudaEventRecord(event_beams[i], stream_beams[i]));
                checkCudaErrors(cudaStreamWaitEvent(stream_beams[0], event_beams[i], 0));
            }
        }
        dose_sum(beams, Phtm, &d_PVCS_total_dose, d_sources, stream_beams[0]);
        beam::calc_FCBB_PVCS_dose_grad(Phtm, &d_element_wise_loss, d_PVCS_total_dose, stream_beams[0]);

        checkCudaErrors(cudaEventRecord(new_event_beams[0], stream_beams[0]));
        for (int i=1; i<stream_beams.size(); i++)
            checkCudaErrors(cudaStreamWaitEvent(stream_beams[i], new_event_beams[0], 0));
        
        for (int i=0; i<beams.size(); i++)
        {
            beams[i].PVCS_dose_backward(Phtm, stream_beams[i]);
            beams[i].BEV_dose_backward(Phtm, FCBB6MeV, stream_beams[i]);
            beams[i].convolveT(FCBB6MeV, stream_beams[i]);
            // beams[i].fluence_map_update(i, )
            
            if (i != 0)
            {
                checkCudaErrors(cudaEventRecord(new_event_beams[i], stream_beams[i]));
                checkCudaErrors(cudaStreamWaitEvent(stream_beams[0], new_event_beams[i], 0));
            }
        }

        reduction(d_element_wise_loss, phantom_size, d_out0, loss, iter, stream_beams.back());
        checkCudaErrors(cudaEventRecord(new_event_beams.back(), stream_beams.back()));
        checkCudaErrors(cudaStreamWaitEvent(stream_beams[0], new_event_beams.back(), 0));
        
        checkCudaErrors(cudaStreamEndCapture(stream_beams[0], &graph_beams));
        checkCudaErrors(cudaGraphInstantiate(&graphExec_beams, graph_beams, NULL, NULL, 0));
        graph_beams_initialized = true;
    }

    cudaStream_t streamForGraph;
    checkCudaErrors(cudaStreamCreate(&streamForGraph));
    for (iter=0; iter<iterations; iter++)
        checkCudaErrors(cudaGraphLaunch(graphExec_beams, streamForGraph));
    cudaStreamSynchronize(streamForGraph);
    checkCudaErrors(cudaMemcpy(h_loss, loss, iterations*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<iterations; i++)
        cout << h_loss[i] << " ";
    cout << endl;
    
    // // Visualize the results
    // optimize_stationary_view_results(beams, Phtm, d_PVCS_total_dose, d_element_wise_loss);
    optimize_stationary_view_backward_results(beams);

    // deconstruct streams and graphs
    checkCudaErrors(cudaGraphExecDestroy(graphExec_beams));
    checkCudaErrors(cudaGraphDestroy(graph_beams));
    for (int i=0; i<beams.size(); i++)
        checkCudaErrors(cudaStreamDestroy(stream_beams[i]));
    checkCudaErrors(cudaStreamDestroy(streamForGraph));
}