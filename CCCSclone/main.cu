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

#include <iostream>
#include <boost/filesystem.hpp>
#include "helper_cuda.h"

namespace fs = boost::filesystem;

old::CONSTANTS* constants;
old::MONO_KERNELS mono_kernels;

int main(int argc, char** argv)
{
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
    std::string full_h5_fname = (fs::path(result_dir) / 
        fs::path("Dose_Coefficients.h5")).string();

    int nrays = constants->nphi * (constants->ntheta / 2);  // convolution directions per beamlet
    int deviceIdx = dev::getarg<int>("deviceIdx");

    std::vector<old::BEAM> beam_arr(constants->beam_count);
    if (old::load_omni_beam_list(beam_arr, constants->beam_count, 0) < 0)
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

    // sanity check
    if (dev::getarg<bool>("debugLog"))
    {
        int flag = old::debugLog(datavols, &mono_kernels, constants, beam_arr);
        if (flag > 0)
            std::cerr << "debug logging failed." << std::endl;
        return 0;
    }

    // prepare results
    // order: beam -> beamlet
    old::RES_LOG result(beam_arr.size());

    // Start calculation
    old::radconvolveTexture(&mono_kernels, constants, beam_arr, nrays, result);

    // write result
    if (writeResults(result))
        return 1;

    old::freeCudaTexture();
    delete constants;
    delete datavols;
    return 0;
}