#include "binary_io.h"
#include "brain_defs.h"
#include "paths.h"
#include "kernel.h"

#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int old::load_density(SHM_DATA* data)
{
    std::string denseFile = (fs::path(Paths::Instance()->temp_dir()) / fs::path("density.raw")).string();
    std::ifstream file(denseFile);
    if (! file)
    {
        std::cerr << "Failed to open density file: " << denseFile << std::endl;
        return false;
    }
    data->density.resize(data->size_data);
    file.read((char*)(data->density.data()), data->size_data*sizeof(float));
    file.close();
    return true;
}

int old::load_data(CONSTANTS* host, SHM_DATA* data)
{
    data->size_data = host->nvoxels();

    if(load_density(data)<0)
    {
        std::cout << "Failed to load density data." << std::endl;
        return false;
    }

    KERNEL kern;
    std::string kern_fname = Paths::Instance()->cum_kernel_file();
    KERNEL::readFromFile(kern, kern_fname);
    data->radial_boundary.resize(kern.nradii);
    memcpy(data->radial_boundary.data(), kern.radial_boundary, kern.nradii*sizeof(float));
    data->kernel_array.resize(kern.ntheta*kern.nradii);
    memcpy(data->kernel_array.data(), kern.total_matrix.data(), kern.ntheta*kern.nradii*sizeof(float));
    host->ntheta = kern.ntheta;
    host->nradii = kern.nradii;

    int numangles = host->nphi * kern.ntheta/2;
    load_convolution_theta_angles(host->conv_theta_deg, numangles);
    load_convolution_phi_angles(host->conv_phi_deg, numangles);

    return 0;
}