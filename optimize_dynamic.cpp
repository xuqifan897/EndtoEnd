#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;
namespace fs = boost::filesystem;

void beam_init(beam* Beam, string& fluence_maps_folder, array<float, 2>& beamAngle, uint beam_idx, phantom& Phtm)
{
    // basic parameter initialization
    (*Beam).zenith = beamAngle[0];
    (*Beam).azimuth = beamAngle[1];
    (*Beam).SAD = get_args<float>("SAD") / 10; // convert mm to cm
    (*Beam).pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;

    vector<int> fluence_map_dimension = get_args<vector<int>>("fluence-map-dimension");
    if (fluence_map_dimension[0]!=E2E::FM_dimension || fluence_map_dimension[1]!=E2E::FM_dimension)
    {
        cout << "Sorry, we only support fluence map \
            dimension of " << E2E::FM_dimension << " at this time" << endl;
        exit;
    }
    (*Beam).fluence_map_dimension = array<int, 2>({fluence_map_dimension[0], fluence_map_dimension[1]});

    vector<int> fluence_map_convolution_radius = \
        get_args<vector<int>>("fluence-map-convolution-radius");
    if (fluence_map_convolution_radius[0]!=E2E::FM_convolution_radius || \
        fluence_map_convolution_radius[1]!=E2E::FM_convolution_radius)
    {
        cout << "Sorry, we only support fluence map \
            convolution radius of " << E2E::FM_convolution_radius << " at this time" << endl;
        exit;
    }

    (*Beam).convolved_fluence_map_dimension = array<int, 2>({ \
        FM_dimension + 2 * FM_convolution_radius, FM_dimension + 2 * FM_convolution_radius});
    (*Beam).extended_fluence_map_dimension = array<int, 2>({ \
        FM_dimension + 4 * FM_convolution_radius, FM_dimension + 4 * FM_convolution_radius});
    
    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    (*Beam).isocenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});

    // log beam parameters
    if (beam_idx == 0)
    {
        cout << "Below are the parameters of the first beam:" << endl;
        cout << "zenith: " << (*Beam).zenith << " rad" << endl;
        cout << "azimuth: " << (*Beam).azimuth << " rad" << endl;
        cout << "SAD: " << (*Beam).SAD << " cm" << endl;
        cout << "pixel size: " << (*Beam).pixel_size << " cm" << endl;
        cout << "fluence map dimension: (" << (*Beam).fluence_map_dimension[0] << \
            ", " << (*Beam).fluence_map_dimension[1] << ")" << endl;
        cout << "convolved fluence map dimension: (" << (*Beam).convolved_fluence_map_dimension[0] << \
            ", " << (*Beam).convolved_fluence_map_dimension[1] << ")" << endl;
        cout << "extended fluence map dimension: (" << (*Beam).extended_fluence_map_dimension[0] << \
            ", " << (*Beam).extended_fluence_map_dimension[1] << ")" << endl;
        cout << "isocenter: (" << (*Beam).isocenter[0] << ", " << (*Beam).isocenter[1] << ", " << \
            (*Beam).isocenter[2] << ")" << endl;
        cout << "sampling range: (" << (*Beam).sampling_range[0] << ", " \
            << (*Beam).sampling_range[1] << ") cm" << endl;
        cout << "number of sampling points: " << (*Beam).sampling_points << endl;
        cout << "\n" << endl;
    }
    
    auto& convolved_fluence_map_dimension = (*Beam).convolved_fluence_map_dimension;
    auto& extended_fluence_map_dimension = (*Beam).extended_fluence_map_dimension;
    checkCudaErrors(cudaMalloc((void**)(&((*Beam).d_convolved_fluence_map)), \
        convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&((*Beam).d_extended_fluence_map)), \
        extended_fluence_map_dimension[0]*extended_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&((*Beam).d_convolved_fluence_map_grad)), \
        convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&((*Beam).d_fluence_grad)), \
        fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));
    
    checkCudaErrors(cudaMalloc((void**)(&((*Beam).d_element_wise_fluence_smoothness_loss)), \
        fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));

    // fluence map initialization
    uint num_elements = (*Beam).extended_fluence_map_dimension[0] * \
            (*Beam).extended_fluence_map_dimension[1];
    float* h_extended_fluence_map = (float*)malloc(num_elements * sizeof(float));
    // nullize all
    for (uint i=0; i<num_elements; i++)
        h_extended_fluence_map[i] = 0;
    if (fluence_maps_folder == string(""))
    {
        for (uint i=2*FM_convolution_radius; i<FM_dimension+2*FM_convolution_radius; i++)
        {
            uint IDX = i * (*Beam).extended_fluence_map_dimension[1];
            for (uint j=2*FM_convolution_radius; j<FM_dimension+2*FM_convolution_radius; j++)
                h_extended_fluence_map[IDX+j] = 1.;
        }
    }
    else
    {
        // string fluence_map_path = fluence_maps_folder + string("/") + format()
        stringstream fluence_map_path_ss;
        fluence_map_path_ss << fluence_maps_folder << "/" << setfill('0') << setw(3) << beam_idx+1 << ".dat";
        string fluence_map_path = fluence_map_path_ss.str();
        cout << "initialize the fluence map from path: " << fluence_map_path << endl;
        ifstream inFile(fluence_map_path);
        if (! inFile.is_open())
        {
            cout << "Could not open the fluence map file: " << fluence_map_path << endl;
            exit;
        }
        float* h_fluence_map = (float*)malloc(fluence_map_dimension[0] \
            * fluence_map_dimension[1] * sizeof(float));
        inFile.read((char*)h_fluence_map, fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float));
        inFile.close();

        // write h_fluence_map to h_extended_fluence_map
        for (uint i=0; i<FM_dimension; i++)
        {
            uint fluence_map_IDX = i * FM_dimension;
            uint extended_fluence_map_IDX = (i + 2 * FM_convolution_radius) * \
                ((*Beam).extended_fluence_map_dimension[1]) + 2 * FM_convolution_radius;
            for (uint j=0; j<FM_dimension; j++)
                h_extended_fluence_map[extended_fluence_map_IDX + j] = h_fluence_map[fluence_map_IDX + j];
        }
        free(h_fluence_map);
    }
    checkCudaErrors(cudaMemcpy((*Beam).d_extended_fluence_map, h_extended_fluence_map, \
        num_elements*sizeof(float), cudaMemcpyHostToDevice));
    
    // // take a look
    // string extended_fluence_map_out_file = get_args<string>("output-folder") + "/" + "extended_fluence_map.dat";
    // ofstream outFile(extended_fluence_map_out_file);
    // if (! outFile.is_open())
    //     cout << "Could not open this file: " << extended_fluence_map_out_file << endl;
    // outFile.write((char*)h_extended_fluence_map, num_elements*sizeof(float));
    // outFile.close();

    free(h_extended_fluence_map);
    (*Beam).FCBBinit(Phtm);
}


int main_for_debug(int argc, char** argv)
{
    // // This is a test function, for fluence map smoothness and its grad
    // if (args_init(argc, argv))
    // {
    //     cerr << "Argument initialization failure." << endl;
    //     exit;
    // }

    // // phantom initialization
    // cout << "\n\n\nPhantom initialization" << endl;
    // phantom Phtm;
    // phantom_init_default(Phtm);
    // Phtm.to_device();
    // Phtm.textureInit();
    
    // // beam initialization
    // beam Beam;
    // string fluence_maps_folder = get_args<string>("fluence-map-init");
    // array<float, 2> beamAngle({2.617995, 0.9424788});
    // beam_init(&Beam, fluence_maps_folder, beamAngle, 0, Phtm);

    // float eta = 0.1;
    // module_test_smoothness_calc(Beam, eta);
    
    module_test_small_reduction();
}


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


void beam_angles_init(vector<array<float, 2>>& beamAngles)
{
    string beam_angle_config_path = get_args<string>("beam-angle-config-path");
    ifstream inFile(beam_angle_config_path);
    if (! inFile.is_open())
    {
        cout << "Could not open the beam angle config file " << beam_angle_config_path << endl;
        exit(EXIT_FAILURE);
    }
    vector<string> lines;
    string line;
    while (getline(inFile, line))
        lines.push_back(line);
    cout << lines.size() << " beams inferred from the beam angle config file." << endl;

    // load to beam beamAngles
    beamAngles.reserve(lines.size());
    for (string& line : lines)
    {
        stringstream ss(line);
        string item;
        vector<string> line_parsed;
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        if (line_parsed.size() != 2)
        {
            cout << "Expected 2 floats per line, but got more or less." << endl;
            exit(EXIT_FAILURE);
        }
        beamAngles.push_back(array<float, 2>({stof(line_parsed[0]), stof(line_parsed[1])}));
    }

    // log beam angles
    cout << "zenith      azimuth" << endl;
    for (uint i=0; i<beamAngles.size(); i++)
        cout << std::fixed << std::setw(7) << std::setprecision(4) << beamAngles[i][0] << "   " << beamAngles[i][1] << endl;
}


void extended_fluence_map_initialization(float** extended_fluence_map, string& fluence_maps_folder, uint idx, float fluence_map_init_value)
{
    uint extended_fluence_map_dimension = FM_dimension + 4 * FM_convolution_radius;
    uint extended_fluence_map_size = extended_fluence_map_dimension * extended_fluence_map_dimension;
    // nullize extended fluence map
    for (uint i=0; i<extended_fluence_map_size; i++)
        (*extended_fluence_map)[i] = 0;
    
    if (fluence_maps_folder == string(""))
    {
        // if fluence_maps_folder is not specified, initialize the fluence map with fluence_map_init_value
        cout << "fluence map initialized to " << fluence_map_init_value << endl;
        for (uint i=0; i<FM_dimension; i++)
        {
            uint IDX = (i + 2 * FM_convolution_radius) * extended_fluence_map_dimension;
            for (uint j=0; j<FM_dimension; j++)
                (*extended_fluence_map)[IDX + j + 2 * FM_convolution_radius] = fluence_map_init_value;
        }
    }
    else
    {
        stringstream fluence_map_path_ss;
        fluence_map_path_ss << fluence_maps_folder << "/" << setfill('0') << setw(3) << idx + 1 << ".dat";
        string fluence_map_path = fluence_map_path_ss.str();
        cout << "fluence map initialized from " << fluence_map_path << endl;

        float* h_fluence_map = (float*)malloc(FM_dimension * FM_dimension * sizeof(float));
        ifstream inFile(fluence_map_path);
        if (! inFile.is_open())
        {
            cout << "Could not open the fluence map path: " << fluence_map_path << endl;
            exit;
        }
        inFile.read((char*)h_fluence_map, FM_dimension * FM_dimension * sizeof(float));
        inFile.close();

        // copy
        for (uint i=0; i<FM_dimension; i++)
        {
            uint fluence_map_idx = i * FM_dimension;
            uint extended_fluence_map_idx = (i + 2 * FM_convolution_radius) * extended_fluence_map_dimension + 2 * FM_convolution_radius;
            for (uint j = 0; j<FM_dimension; j++)
                (*extended_fluence_map)[extended_fluence_map_idx + j] = h_fluence_map[fluence_map_idx + j];
        }
        free(h_fluence_map);
    }
}


void beams_init_optimize_stationary(vector<beam>& beams, phantom& Phtm)
{
    // static initialization
    beam::FCBBStaticInit(Phtm);
    string fluence_maps_folder = get_args<string>("fluence-map-init");
    
    // initialize beam angles
    vector<array<float, 2>> beamAngles;
    beam_angles_init(beamAngles);

    // fluence_map_init_value is used to initialize the extended fluence map when fluence-map-init is not specified
    float fluence_map_init_value = get_args<float>("ideal-dose") / (beamAngles.size() * get_args<float>("isocenter-dose"));
    if (fluence_maps_folder == string(""))
        cout << "\nfluence maps initialization folder not found. Fluence maps are initialized to " << fluence_map_init_value << "." << endl;
    else
        cout << "\nfluence maps initialization folder found. Fluence maps would be initialized from the binary files in it." << endl;

    vector<int> fluence_map_dimension = get_args<vector<int>>("fluence-map-dimension");
    if (fluence_map_dimension[0]!=E2E::FM_dimension || fluence_map_dimension[1]!=E2E::FM_dimension)
    {
        cout << "Sorry, we only support fluence map \
            dimension of " << E2E::FM_dimension << " at this time" << endl;
        exit;
    }

    vector<int> fluence_map_convolution_radius = \
        get_args<vector<int>>("fluence-map-convolution-radius");
    if (fluence_map_convolution_radius[0]!=E2E::FM_convolution_radius || \
        fluence_map_convolution_radius[1]!=E2E::FM_convolution_radius)
    {
        cout << "Sorry, we only support fluence map \
            convolution radius of " << E2E::FM_convolution_radius << " at this time" << endl;
        exit;
    }

    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    // convert mm to cm
    isocenter[0] /= 10;
    isocenter[1] /= 10;
    isocenter[2] /= 10;

    beams.reserve(beamAngles.size());
    uint extended_fluence_map_size = (FM_dimension + 4 * FM_convolution_radius) * \
        (FM_dimension + 4 * FM_convolution_radius);
    
    // the memory to store extended fluence map
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size * sizeof(float));

    for (uint i=0; i<beamAngles.size(); i++)
    {
        beams.push_back(beam());
        beam& new_beam = beams.back();
        new_beam.zenith = beamAngles[i][0];
        new_beam.azimuth = beamAngles[i][1];
        new_beam.SAD = get_args<float>("SAD") / 10; // convert mm to cm
        new_beam.pixel_size = get_args<vector<float>>("fluence-map-pixel-size")[0] / 10;
        new_beam.fluence_map_dimension = array<int, 2>({FM_dimension, FM_dimension});
        new_beam.convolved_fluence_map_dimension = array<int, 2>({ \
            FM_dimension + 2 * FM_convolution_radius, FM_dimension + 2 * FM_convolution_radius});
        new_beam.extended_fluence_map_dimension = array<int, 2>({ \
            FM_dimension + 4 * FM_convolution_radius, FM_dimension + 4 * FM_convolution_radius});
        new_beam.isocenter = array<float, 3>({isocenter[0], isocenter[1], isocenter[2]});

        // log beam parameters
        cout << "Beam " << i + 1 << " parameters:" << endl;
        cout << "zenith: " << new_beam.zenith << " rad" << endl;
        cout << "azimuth: " << new_beam.azimuth << " rad" << endl;
        cout << "SAD: " << new_beam.SAD << " cm" << endl;
        cout << "pixel size: " << new_beam.pixel_size << " cm" << endl;
        cout << "fluence map dimension: (" << new_beam.fluence_map_dimension[0] << \
            ", " << new_beam.fluence_map_dimension[1] << ")" << endl;
        cout << "convolved fluence map dimension: (" << new_beam.convolved_fluence_map_dimension[0] << \
            ", " << new_beam.convolved_fluence_map_dimension[1] << ")" << endl;
        cout << "extended fluence map dimension: (" << new_beam.extended_fluence_map_dimension[0] << \
            ", " << new_beam.extended_fluence_map_dimension[1] << ")" << endl;
        cout << "isocenter: (" << new_beam.isocenter[0] << ", " << new_beam.isocenter[1] << ", " << \
            new_beam.isocenter[2] << ")" << endl;
        cout << "sampling range: (" << new_beam.sampling_range[0] << ", " \
            << new_beam.sampling_range[1] << ") cm" << endl;
        cout << "number of sampling points: " << new_beam.sampling_points << endl;

        // cuda memory allocation
        auto& convolved_fluence_map_dimension = new_beam.convolved_fluence_map_dimension;
        auto& extended_fluence_map_dimension = new_beam.extended_fluence_map_dimension;
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_convolved_fluence_map)), \
            convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_extended_fluence_map)), \
            extended_fluence_map_dimension[0]*extended_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_convolved_fluence_map_grad)), \
            convolved_fluence_map_dimension[0]*convolved_fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_fluence_grad)), \
            fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));
        checkCudaErrors(cudaMalloc((void**)(&(new_beam.d_element_wise_fluence_smoothness_loss)), \
            fluence_map_dimension[0]*fluence_map_dimension[1]*sizeof(float)));
        
        // extended fluence map initialization
        extended_fluence_map_initialization(&h_extended_fluence_map, fluence_maps_folder, i, fluence_map_init_value);
        checkCudaErrors(cudaMemcpy(new_beam.d_extended_fluence_map, h_extended_fluence_map, \
            extended_fluence_map_size*sizeof(float), cudaMemcpyHostToDevice));
        
        new_beam.FCBBinit(Phtm);
        cout << "\n" << endl;
    }
    free(h_extended_fluence_map);
}


void optimize(vector<beam>& beams, phantom& Phtm, FCBBkernel* kernel, float** h_loss, float** h_smoothness_loss)
{
    int iterations = get_args<int>("iterations");
    float step_size = get_args<float>("step-size");
    float eta = get_args<float>("eta");

    // initialization of necessary cuda arrays
    uint phantom_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* d_PVCS_total_dose;
    float* d_element_wise_loss;
    checkCudaErrors(cudaMalloc((void**)(&d_PVCS_total_dose), phantom_size*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)(&d_element_wise_loss), phantom_size*sizeof(float)));
    
    float* d_out0 = nullptr; // intermediate result for reduction
    float* loss = nullptr; // the final phantom dose loss array
    checkCudaErrors(cudaMalloc((void**)&d_out0, REDUCTION_BLOCK_SIZE*sizeof(float)));
    checkCudaErrors(cudaMalloc((void**)&loss, iterations*sizeof(float)));

    float* smoothness_loss = nullptr;
    checkCudaErrors(cudaMalloc((void**)(&smoothness_loss), iterations*beams.size()*sizeof(float)));

    // for fluence map update
    float** d_squared_grad = (float**)malloc(beams.size() * sizeof(float*));
    for (uint i=0; i<beams.size(); i++)
        checkCudaErrors(cudaMalloc(d_squared_grad+i, FM_dimension*FM_dimension*sizeof(float)));
    float* d_norm_final;
    checkCudaErrors(cudaMalloc(&d_norm_final, beams.size()*sizeof(float)));

    // d_sources and h_sources are the arrays to contain the d_FCBB_PVCS_dose
    float** h_sources = (float**)malloc(beams.size()*sizeof(float*));
    float** d_sources = nullptr;
    for (int i=0; i<beams.size(); i++)
    {
        if (beams[i].d_FCBB_PVCS_dose == nullptr)
        {
            cout << "The " << i << "th beam has not called FCBBinit()" << endl;
            exit(EXIT_FAILURE);
        }
        h_sources[i] = beams[i].d_FCBB_PVCS_dose;
    }
    checkCudaErrors(cudaMalloc((void***)&d_sources, beams.size()*sizeof(float*)));
    checkCudaErrors(cudaMemcpy(d_sources, h_sources, beams.size()*sizeof(float*), \
        cudaMemcpyHostToDevice));

    // h_loss and h_smoothness_loss initialization
    if (*h_loss == nullptr)
        *h_loss = (float*)malloc(iterations*sizeof(float));
    if (*h_smoothness_loss == nullptr)
        *h_smoothness_loss = (float*)malloc(iterations*beams.size()*sizeof(float));
    
    // allocate perturbation loss

    // cuda event for time measurement
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    for (int iter = 0; iter<iterations; iter++)
    {
        // forward pass: calculate each beam's dose
        cudaEventRecord(start);
        for (int i=0; i<beams.size(); i++)
        {
            beams[i].convolve(kernel);
            beams[i].BEV_dose_forward(Phtm, FCBB6MeV);
            beams[i].PVCS_dose_forward(Phtm);
        }

        dose_sum(beams, Phtm, &d_PVCS_total_dose, d_sources);
        beam::calc_FCBB_PVCS_dose_grad(Phtm, &d_element_wise_loss, d_PVCS_total_dose);
        reduction(d_element_wise_loss, phantom_size, d_out0, loss, iter);

        // update fluence map
        for (int i=0; i<beams.size(); i++)
        {
            beams[i].PVCS_dose_backward(Phtm);
            beams[i].BEV_dose_backward(Phtm, FCBB6MeV);
            beams[i].convolveT(FCBB6MeV);

            // smoothness update
            beams[i].smoothness_calc(eta);
            uint loss_idx = iter * beams.size() + i;
            reduction_small(beams[i].d_element_wise_fluence_smoothness_loss, \
                FM_dimension*FM_dimension, smoothness_loss, loss_idx);

            beams[i].fluence_map_update(i, d_norm_final, d_squared_grad[i], step_size);
        }
        cudaEventRecord(stop);
        // log
        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        checkCudaErrors(cudaMemcpy(*h_loss, loss, iterations*sizeof(float), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(*h_smoothness_loss, smoothness_loss, \
            iterations*beams.size()*sizeof(float), cudaMemcpyDeviceToHost));
        float total_smoothness_loss = 0;
        for (int i=0; i<beams.size(); i++)
        {
            int idx = iter * beams.size() + i;
            total_smoothness_loss += (*h_smoothness_loss)[idx];
        }
        cout << "Iteration: " << iter + 1 << ", step size: " << step_size << \
            ", execution time: " << milliseconds << " ms, dose loss: " << (*h_loss)[iter] << \
            ", smoothness loss: " << total_smoothness_loss << endl;
    }

    // clean up
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    checkCudaErrors(cudaFree(d_PVCS_total_dose));
    checkCudaErrors(cudaFree(d_element_wise_loss));
    checkCudaErrors(cudaFree(d_out0));
    checkCudaErrors(cudaFree(loss));
    checkCudaErrors(cudaFree(smoothness_loss));
    for (int i=0; i<beams.size(); i++)
        checkCudaErrors(cudaFree(d_squared_grad[i]));
    free(d_squared_grad);
    checkCudaErrors(cudaFree(d_norm_final));
    checkCudaErrors(cudaFree(d_sources));
    free(h_sources);
}


void extended_to_fluence(float* h_fluence_map, float* h_extended_fluence_map)
{
    uint extended_fluence_map_dimension = FM_dimension + 4 * FM_convolution_radius;
    for (uint i=0; i<FM_dimension; i++)
    {
        uint fluence_map_idx = i * FM_dimension;
        uint extended_fluence_map_idx = (i + 2 * FM_convolution_radius) * \
            extended_fluence_map_dimension + 2 * FM_convolution_radius;
        for (uint j=0; j<FM_dimension; j++)
            h_fluence_map[fluence_map_idx + j] = h_extended_fluence_map[extended_fluence_map_idx + j];
    }
}


void log_out(vector<beam>& beams, phantom& Phtm, FCBBkernel* kernel, float* h_loss, float* h_smoothness_loss)
{
    string output_folder = get_args<string>("output-folder");
    fs::path output_folder_fs(output_folder);
    if (! fs::is_directory(output_folder_fs))
        fs::create_directory(output_folder_fs);
    
    string fluence_map_folder = output_folder + "/" + "fluence_maps";
    fs::path fluence_map_folder_fs(fluence_map_folder);
    if (! fs::is_directory(fluence_map_folder_fs))
        fs::create_directory(fluence_map_folder_fs);
    
    // output individual fluence maps
    uint extended_fluence_map_dimension = FM_dimension + 4 * FM_convolution_radius;
    uint extended_fluence_map_size = extended_fluence_map_dimension * \
        extended_fluence_map_dimension;
    float* h_extended_fluence_map = (float*)malloc(extended_fluence_map_size * sizeof(float));
    float* h_fluence_map = (float*)malloc(FM_dimension * FM_dimension * sizeof(float));

    for (int idx=0; idx<beams.size(); idx++)
    {
        beam& this_beam = beams[idx];
        checkCudaErrors(cudaMemcpy(h_extended_fluence_map, this_beam.d_extended_fluence_map, \
            extended_fluence_map_size*sizeof(float), cudaMemcpyDeviceToHost));
        extended_to_fluence(h_fluence_map, h_extended_fluence_map);
        
        stringstream outputPath_ss;
        outputPath_ss << fluence_map_folder << "/" << setfill('0') << setw(3) << idx + 1 << ".dat";
        string outputPath = outputPath_ss.str();
        ofstream outFile(outputPath);
        if (! outFile.is_open())
        {
            cout << "Could not open the output fluence map file: " << outputPath << endl;
            exit(EXIT_FAILURE);
        }
        outFile.write((char*)h_fluence_map, FM_dimension*FM_dimension*sizeof(float));
        outFile.close();
        cout << "Fluence map " << idx + 1 << " writtent to " << outputPath << endl;
    }

    // output total dose
    uint dose_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_PVCS_total_dose = (float*)malloc(dose_size*sizeof(float));
    float* d_PVCS_total_dose;
    checkCudaErrors(cudaMalloc((void**)(&d_PVCS_total_dose), dose_size*sizeof(float)));
    for (int idx=0; idx<beams.size(); idx++)
    {
        beam& this_beam = beams[idx];
        this_beam.convolve(kernel);
        this_beam.BEV_dose_forward(Phtm, kernel);
        this_beam.PVCS_dose_forward(Phtm);
    }
    float** h_sources = (float**)malloc(beams.size()*sizeof(float*));
    for (int idx=0; idx<beams.size(); idx++)
        h_sources[idx] = beams[idx].d_FCBB_PVCS_dose;
    float ** d_sources;
    checkCudaErrors(cudaMalloc((void***)(&d_sources), beams.size()*sizeof(float*)));
    checkCudaErrors(cudaMemcpy(d_sources, h_sources, beams.size()*sizeof(float*), cudaMemcpyHostToDevice));
    dose_sum(beams, Phtm, &d_PVCS_total_dose, d_sources);
    checkCudaErrors(cudaMemcpy(h_PVCS_total_dose, d_PVCS_total_dose, dose_size*sizeof(float), cudaMemcpyDeviceToHost));
    string outputPath = output_folder + "/totalDose.dat";
    ofstream outFile(outputPath);
    if (! outFile.is_open())
    {
        cout << "Could not open total dose file: " << outputPath << endl;
        exit(EXIT_FAILURE);
    }
    outFile.write((char*)h_PVCS_total_dose, dose_size*sizeof(float));
    outFile.close();
    cout << "Total dose written to " << outputPath << endl;

    // output loss function
    string lossOutputPath = output_folder + "/DoseLoss.dat";
    outFile.open(lossOutputPath);
    if (! outFile.is_open())
    {
        cout << "Could not open dose loss output file: " << lossOutputPath << endl;
        exit(EXIT_FAILURE);
    }
    outFile.write((char*)h_loss, get_args<int>("iterations")*sizeof(float));
    outFile.close();
    cout << "Dose loss written to " << lossOutputPath << endl;

    lossOutputPath = output_folder + "/SmoothnessLoss.dat";
    outFile.open(lossOutputPath);
    if (! outFile.is_open())
    {
        cout << "Could not open smoothness loss output file: " << lossOutputPath << endl;
        exit(EXIT_FAILURE);
    }
    outFile.write((char*)h_smoothness_loss, get_args<int>("iterations")*beams.size()*sizeof(float));
    outFile.close();
    cout << "Smoothness loss written to " << lossOutputPath << endl;

    // clean up
    free(h_fluence_map);
    free(h_extended_fluence_map);
    free(h_PVCS_total_dose);
    checkCudaErrors(cudaFree(d_PVCS_total_dose));
    free(h_sources);
    checkCudaErrors(cudaFree(d_sources));
}


int main(int argc, char** argv)
{
    // This functino optimizes fluence map alone, without optimizing beam angles.
    // Fluence map smoothness term is included
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
    vector<beam> beams;
    beams_init_optimize_stationary(beams, Phtm);

    // begin optimize
    float* h_loss = nullptr; // Dose losses, which will be allocated inside optimize function
    float* h_smoothness_loss = nullptr; // Fluence map smoothness loss, which will be allocated inside the optimize function
    optimize(beams, Phtm, kernel, &h_loss, &h_smoothness_loss);

    // log out
    log_out(beams, Phtm, kernel, h_loss, h_smoothness_loss);

    // clean up
    free(h_loss);
    free(h_smoothness_loss);
}