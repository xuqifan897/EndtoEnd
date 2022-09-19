#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;


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

void beam_dose_calculation(vector<array<float, 2>>& beamAngles, phantom& Phtm, FCBBkernel* kernel)
{
    // static initialization
    beam::FCBBStaticInit(Phtm);
    string fluence_maps_folder = get_args<string>("fluence-map-init");
    if (fluence_maps_folder == string(""))
        cout << "fluence maps initialization folder not found. Fluence maps are initialized to one.\n" << endl;
    beam* Beam;

    uint dose_size = Phtm.dimension[0] * Phtm.dimension[1] * Phtm.pitch;
    float* h_dose = (float*)malloc(dose_size*sizeof(float));

    for (uint i=0; i<beamAngles.size(); i++)
    {
        // beam initialization
        Beam = new beam;
        array<float, 2>& beamAngle = beamAngles[i];
        beam_init(Beam, fluence_maps_folder, beamAngle, i, Phtm);

        // dose calculation
        (*Beam).convolve(kernel);
        (*Beam).BEV_dose_forward(Phtm, kernel);
        (*Beam).PVCS_dose_forward(Phtm);

        // 
        checkCudaErrors(cudaMemcpy(h_dose, (*Beam).d_FCBB_PVCS_dose, dose_size*sizeof(float), cudaMemcpyDeviceToHost));
        stringstream doseOutputPath_ss;
        doseOutputPath_ss << get_args<string>("output-folder") << "/dose" << setfill('0') << setw(3) << i+1 << ".dat";
        string doseOutputPath = doseOutputPath_ss.str();
        ofstream outFile(doseOutputPath);
        if (! outFile.is_open())
        {
            cout << "Could not open output file: " << doseOutputPath << endl;
            exit;
        }
        outFile.write((char*)h_dose, dose_size*sizeof(float));
        outFile.close();
        cout << "result dose written to " << doseOutputPath << endl;

        // cleanup
        delete Beam;
    }
    free(h_dose);
    beam::FCBBStaticDecon();
}

int main(int argc, char** argv)
{
    // this function produces the dose for each beam in a beam list file
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

    // initialize beam angles
    vector<array<float, 2>> beamAngles;
    beam_angles_init(beamAngles);

    // beam dose calculation
    cout << "\n\n\nbeam dose calculation" << endl;
    beam_dose_calculation(beamAngles, Phtm, kernel);
}