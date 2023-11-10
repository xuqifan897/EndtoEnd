#include "kernel.h"
#include "paths.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "H5Cpp.h"

// fillers for these entries in kernel structure
#define UNCERT 0.0
#define MEAN_RADIUS 0.0
#define MEAN_ANGLE 0.0

int old::KERNEL::readFromFile(KERNEL& kern, std::string fname, bool verbose)
{
    H5::Exception::dontPrint();
    try
    {
        H5::H5File h5file(fname, H5F_ACC_RDONLY);
        H5::Group rootgroup = h5file.openGroup("/");
        if (!KERNEL::_readFromHDF5(kern, rootgroup))
        {
            if (verbose)
                std::cout << "Failed to read KERNEL from " << fname << std::endl;
            return false;
        }
    }
    catch (H5::FileIException &file_exists_error)
    {
        if (verbose)
            std::cout << "Failed to read KERNEL from " << fname << std::endl;
        return false;
    }
    return true;
}

int old::KERNEL::_readFromHDF5(KERNEL& kern, H5::Group& h5group)
{
    {
        //read kernel_file string
        auto att = h5group.openAttribute("kernel_file");
        H5std_string buff("");
        att.read(att.getDataType(), buff);
        kern.kernel_file = buff;
    }
    {
        // read angles, radii
        auto att = h5group.openAttribute("nradii");
        att.read(H5::PredType::NATIVE_INT, &kern.nradii);
        att = h5group.openAttribute("ntheta");
        att.read(H5::PredType::NATIVE_INT, &kern.ntheta);

        att = h5group.openAttribute("radial_boundaries_cm");
        att.read(H5::PredType::NATIVE_FLOAT, &kern.radial_boundary);

        att = h5group.openAttribute("angular_boundaries_deg");
        att.read(H5::PredType::NATIVE_FLOAT, &kern.angular_boundary);
    }
    {
        // read matrix
        auto dset = h5group.openDataSet("weights");
        auto file_space = dset.getSpace();
        hsize_t dims[3] {};
        file_space.getSimpleExtentDims(dims);
        if (dims[0] != N_KERNEL_CATEGORIES || dims[1] != kern.ntheta || dims[2] != kern.nradii)
        {
            throw std::runtime_error("dataset dimensions do not agree with kernel dimensions");
        }
        hsize_t subcount[] = {1, uint(kern.ntheta), uint(kern.nradii)};
        auto mem_space = H5::DataSpace(3, subcount);
        for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
            hsize_t offset[] = {uint(cc), 0, 0};
            kern.matrix[cc] = std::vector<float>(kern.nradii * kern.ntheta);
            file_space.selectHyperslab(H5S_SELECT_SET, subcount, offset);
            dset.read(kern.matrix[cc].data(), H5::PredType::NATIVE_FLOAT, mem_space, file_space);
        }
    }
    {
        //read total matrix
        auto dset = h5group.openDataSet("total_weights");
        auto file_space = dset.getSpace();
        hsize_t dims[2] {};
        file_space.getSimpleExtentDims(dims);
        if (dims[0] != kern.ntheta || dims[1] != kern.nradii)
        {
            throw std::runtime_error("dataset dimensions do not agree with kernel dimensions");
        }
        kern.total_matrix = std::vector<float>(kern.nradii * kern.ntheta);
        dset.read(kern.total_matrix.data(), H5::PredType::NATIVE_FLOAT, file_space, file_space);
    }
    return true;
}

int _load_conv_angles(std::vector<float>& result, int size, std::string fname)
{
    std::ifstream File(fname);
    if (! File)
    {
        std::cout << "Failed to open convolution angles file " << fname << std::endl;
        return false;
    }
    result.resize(size);
    File.read((char*)(result.data()), size*sizeof(float));
    File.close();
    return true;
}

int old::load_convolution_theta_angles(std::vector<float>& conv_theta_deg, int size)
{
    bool result = _load_conv_angles(conv_theta_deg, size, Paths::Instance()->conv_theta_file());
    return result;
}

int old::load_convolution_phi_angles(std::vector<float>& conv_phi_deg, int size)
{
    bool result = _load_conv_angles(conv_phi_deg, size, Paths::Instance()->conv_phi_file());
    return result;
}