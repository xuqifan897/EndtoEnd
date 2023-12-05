/* this code is to benchmark the runtime of the existing dose calculation algorithm */
#include "argparse.h"
#include "beam.h"
#include "brain_defs.h"
#include "calc_header.h"
#include "binary_io.h"
#include "paths.h"
#include "kernel.h"
#include "spectra.h"
#include "cudaInit.h"
#include "cudaDoseCalc.h"
#include "debugLog.h"

#include "benchmark_host.cuh"

#include <iostream>
#include <boost/filesystem.hpp>
#include "helper_cuda.h"

namespace fs = boost::filesystem;

old::CONSTANTS* constants;
old::MONO_KERNELS mono_kernels;

int main(int argc, char** argv) {
    if (dev::argparse(argc, argv))
        return 0;
    
    // load pre computed data
    constants = new old::CONSTANTS;
    old::constantsInit(constants);

    old::SHM_DATA* datavols = new old::SHM_DATA;
    old::load_data(constants, datavols);

    // io purpose
    std::string result_dir = old::Paths::Instance()->result_dir();
    if (! fs::is_directory(result_dir))
        fs::create_directories(result_dir);
    
    int nrays = constants->nphi * (constants->ntheta / 2);
    int deviceIdx = dev::getarg<int>("deviceIdx");

    std::vector<old::BEAM> beam_arr;
    if (dev::beamsInit(beam_arr) < 0)
    {
        std::cerr << "Failed to load beam list." << std::endl;
        return -1;
    }

    // load spectrum files
    mono_kernels.spectrum_file = std::string(constants->beam_spec);

    // monoenergetic kernel file
    // spectrum data
    if (read_spectrum_file(&mono_kernels, 1) != 1)
    {
        std::cerr << "Failed reading spectrum file!\n" << std::endl;
        return -1;
    }

    // Store copies of constant problem data and texture / surface references on GPU
    std::cout << "Initializing memory on GPU: " << deviceIdx << std::endl;
    checkCudaErrors(cudaSetDevice(deviceIdx));
    old::initCudaConstandTex(datavols, &mono_kernels, constants, nrays);

    dev::radconvTextureBench(&mono_kernels, constants, beam_arr, nrays);
}