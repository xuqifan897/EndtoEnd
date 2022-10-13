#include <boost/program_options.hpp>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include "args.h"

namespace po = boost::program_options;
using namespace std;
po::variables_map* E2E::args=nullptr;
uint E2E::pitch_module=8;

void log_arguments(po::variables_map& vm);
void iterate_args_str();

int E2E::args_init(int argc, char** argv)
{   
    try
    {   
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("phantom-dimension", po::value<vector<int>>()->multitoken(), "the dimension of the phantom, of size 3")
            ("voxel-size", po::value<float>(), "the voxel size of the phantom in mm. Each voxel dimension is equal.")
            ("phantom-isocenter", po::value<vector<float>>()->multitoken(), "the isocenter coordinates in mm. In x, y, z order")
            ("phantom-path", po::value<string>(), "the file path to the phantom. The file is a binary file, \
            //     with voxels are ranged in order x, y, z, the same as below. Example is generated from numpy.array.totile() method")
            ("PTV-weight-path", po::value<string>(), "PTV weight file path, refer to README for detail")
            ("PTV-target-path", po::value<string>(), "PTV target file path, refer to README for detail")
            ("OAR-weight-path", po::value<string>(), "OAR weight file path, refer to README for detail")
            ("OAR-target-path", po::value<string>(), "OAR target file path, refer to README for detail")
            ("fluence-map-init", po::value<string>()->default_value(""), "The folder containing fluence map initializations. \
                If not specified, the fluence map is initialized to all ones")
            ("beam-energy", po::value<float>()->default_value(6.), "the energy of the beam in MeV. At this point we only \
                support monoenergetic beam. The beam energy should be chosen from 4, 6, 10, 15, 24")
            ("SAD", po::value<float>()->default_value(1000.), "Source-to-axis distance, in mm")
            ("number-of-beams", po::value<int>()->default_value(30), "Number of beams used in the optimization")
            ("fluence-map-dimension", po::value<vector<int>>()->multitoken(), "Fluence map dimension \
                at reference plane, of size 2")
            ("fluence-map-convolution-radius", po::value<vector<int>>()->multitoken(), "The radius of fluence map convolution kernel, in # pixel")
            ("fluence-map-sampling-range", po::value<vector<float>>()->multitoken(), "The sampling range along a ray in mm, of size 2")
            ("fluence-map-sampling-points", po::value<int>()->default_value(512), "the number of sampling points along a ray")
            ("fluence-map-pixel-size", po::value<vector<float>>()->multitoken(), "the pixel size of fluence map")
            // ("fluence-map-output-path", po::value<string>(), "The output path of the fluence maps, in order n, x, y, where \
            //     n is the number of beams")
            ("output-folder", po::value<string>(), "Output folder path")
            ("zenith-range", po::value<vector<float>>()->multitoken(), "The zenith range of the beam to avoid physical collision")
            // ("dose-path", po::value<string>(), "The path to which the dose is stored")
            ("spectrum-path", po::value<string>(), "The path to spectrum data")
            ("ATheta-path", po::value<string>(), "The path to ATheta file")
            ("atheta-path", po::value<string>(), "The path to atheta file")
            ("BTheta-path", po::value<string>(), "The path to BTheta file")
            ("btheta-path", po::value<string>(), "The path to btheta file")
            ("pencil-path", po::value<string>(), "The path to the pencil beam lateral kernel file")
            ("depthDose-path", po::value<string>(), "The path to the depth dose table file")
            ("beam-angle-config-path", po::value<string>(), "The path to the beam angle configuration file")
            ("iterations", po::value<int>()->default_value(100), "The number of iterations in the one stage. Three stages in total, with step size scaling.")
            ("step-size", po::value<float>()->default_value(5e-3), "fluence map update step size for the first stage.")
            ("zenith", po::value<float>()->default_value(PI/2), "The zenith angle of the beam, used for single beam dose calculation")
            ("azimuth", po::value<float>()->default_value(0), "The azimuth angle of the beam, used for single beam dose calculation")
            ("ideal-dose", po::value<float>()->default_value(20.), "The expected PTV dose, used for initialization")
            ("isocenter-dose", po::value<float>()->default_value(3.54), "Emperical value of the isocenter dose when the fluence map is \
                set to 1. This value is used to initialize the fluence map if fluence-map-init is not initialized")
            ("eta", po::value<float>()->default_value(1e3), "eta, the weight for fluence map smoothness.")
            ("step-size-angular", po::value<float>()->default_value(0.03), "Step size for zenith and azimuth angle updates")
            ("step-size-angular-max", po::value<float>()->default_value(0.6), "maximum step size for zenith and azimuth angle updates, which is linearly ramped down")
            ("step-size-angular-min", po::value<float>()->default_value(0.06), "minimum step size for zenith and azimuth angle updates, which is linearly ramped down")
            ("temperature", po::value<float>()->default_value(1e5), "temperature for simulated annealing")
            ("temperature-max", po::value<float>()->default_value(4e5), "maximum temperature for simulated annealing. The actual temperature is linearly ramped down")
            ("temperature-min", po::value<float>()->default_value(1e5), "minimum temperature for simulated annealing. The actual temperature is linearly ramped down")
            ("random-seed", po::value<uint>()->default_value(10086), "random number seed for annealing initialization")
        ;

        E2E::args = new po::variables_map();
        // E2E::args = (po::variables_map*)malloc(sizeof(po::variables_map));
        // (*E2E::args) = po::variables_map();
        po::store(po::parse_command_line(argc, argv, desc), *E2E::args);
        po::notify(*E2E::args);

        if ((*E2E::args).count("help"))
        {
            cout << desc << endl;
            return 1;
        }

        log_arguments(*E2E::args);
        if (spectrum_init()) return 1;
        if (CCCSkernel_init()) return 1;
        if (FCBBkernel_init()) return 1;
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }    
    return 0;
}

void log_arguments(po::variables_map& vm)
{
    int num_chars = 60;
    for(const auto& it : vm)
    {   
        string output = it.first;
        string output_{};
        try
        {
            output_ = vm[output].as<string>();
            output = output + string(10, '.') + output_;
            cout << output << endl;
            continue;
        }
        catch(const std::exception& e){}

        try
        {
            int v = vm[output].as<int>();
            output_ = to_string(v);
            goto HERE;
        }
        catch(const std::exception& e){}

        try
        {
            float v = vm[output].as<float>();
            output_ = to_string(v);
            goto HERE;
        }
        catch(const std::exception& e){}

        try
        {
            vector<int> v = vm[output].as<vector<int>>();
            if (v.size() == 3)
            {
                output_ = "(" + to_string(v[0]) + ", " + to_string(v[1]) + ", " + to_string(v[2]) + ")";
            }
            else if(v.size() == 2)
            {
                output_ = "(" + to_string(v[0]) + ", " + to_string(v[1]) + ")";
            }
            goto HERE;
        }
        catch(const std::exception& e){}

        try
        {
            vector<float> v = vm[output].as<vector<float>>();
            if (v.size() == 3)
            {
                output_ = "(" + to_string(v[0]) + ", " + to_string(v[1]) + ", " + to_string(v[2]) + ")";
            }
            else if(v.size() == 2)
            {
                output_ = "(" + to_string(v[0]) + ", " + to_string(v[1]) + ")";
            }
            goto HERE;
        }
        catch(const std::exception& e){}
        
        HERE: int num_dots = num_chars - output.size() - output_.size();
        cout << output + string(max(num_dots, 10), '.') + output_ << endl;
    }
}

void E2E::deviceProperty()
{
    int nDevices;
    int deviceSupportsMemoryPools;
    cudaGetDeviceCount(&nDevices);
    for (int i=0; i<nDevices; i++)
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        cout << "Global memory: " << prop.totalGlobalMem << endl;
        cout << "Shared memory per block: " << prop.sharedMemPerBlock << endl;
        cout << "Shared memory per multiprocessor: " << prop.sharedMemPerMultiprocessor << endl;
        cout << "Warp size: " << prop.warpSize << endl;
        cout << "Max threads per block: " << prop.maxThreadsPerBlock << endl;
        cout << "Max threads per multiprocessor: " << prop.maxThreadsPerMultiProcessor << endl;
        cout << "Number of multiprocessors: " << prop.multiProcessorCount << endl;
        cout << "Max thread dimension: (" << prop.maxThreadsDim[0] << ", " \
            << prop.maxThreadsDim[1] << ", " << prop.maxThreadsDim[2] << ")" << endl;
        
        // checkCudaErrors(cudaDeviceGetArribute(&deviceSupportsMemoryPools, \
        //     cudaDevAttrMemoryPoolsSupported, i));
    }
}