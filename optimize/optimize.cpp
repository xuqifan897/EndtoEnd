#include <vector>
#include <array>
#include <string>
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

void optimize_stationary_view_results(vector<beam>& beams, phantom& Phtm)
{
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
}

void E2E::optimize_stationary(vector<beam>& beams, phantom& Phtm)
{
    // initialize the extended fluence map of beams
    // the valid fluence map zone is set to 1
    fluence_map_init(beams, Phtm);

    // intermediate data initialization
    float* d_total_dose = nullptr;

    // initialize streams and graphs
    vector<cudaStream_t> stream_beams(beams.size());
    vector<cudaEvent_t> event_beams(beams.size());
    cudaGraph_t graph_beams;
    cudaGraphExec_t graphExec_beams;
    bool graph_beams_initialized = false;
    for (int i=0; i<beams.size(); i++)
        checkCudaErrors(cudaStreamCreate(&(stream_beams[i])));
    for (int i=0; i<beams.size(); i++)
        checkCudaErrors(cudaEventCreate(&(event_beams[i])));
    
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
        
        checkCudaErrors(cudaStreamEndCapture(stream_beams[0], &graph_beams));
        checkCudaErrors(cudaGraphInstantiate(&graphExec_beams, graph_beams, NULL, NULL, 0));
        graph_beams_initialized = true;
    }

    cudaStream_t streamForGraph;
    checkCudaErrors(cudaStreamCreate(&streamForGraph));
    // int iterations = get_args<float>("iterations");
    // for (int i=0; i<iterations; i++)
    // {
    //     checkCudaErrors(cudaGraphLaunch(graphExec_beams, streamForGraph));
    // }
    
    // for debug purposes, here we only launch the graph for 1 iteration, 
    // and visualize the results
    checkCudaErrors(cudaGraphLaunch(graphExec_beams, streamForGraph));
    cudaStreamSynchronize(streamForGraph);
    optimize_stationary_view_results(beams, Phtm);

    // deconstruct streams and graphs
    checkCudaErrors(cudaGraphExecDestroy(graphExec_beams));
    checkCudaErrors(cudaGraphDestroy(graph_beams));
    for (int i=0; i<beams.size(); i++)
        checkCudaErrors(cudaStreamDestroy(stream_beams[i]));
    checkCudaErrors(cudaStreamDestroy(streamForGraph));
}