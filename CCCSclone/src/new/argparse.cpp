#include "argparse.h"
#include "utils.h"
#include <boost/program_options.hpp>
#include <string>
#include <array>
#include <iostream>
namespace po = boost::program_options;

po::variables_map dev::vm;
std::vector<float> dev::dicomVolumeStartCoords;

bool dev::argparse(int argc, char** argv)
{
    po::options_description desc("The argument list for the new CCCS algorithm validation");
    desc.add_options()
        ("help", "Produce help messages. This program requires the following files in the data folder: \n"
            "1) convolution_phi_angles.raw, which contains the convolution phi angles for each convolution direction.\n"
            "2) convolution_theta_angles.raw, which contains the convolution theta angles for each convolution direction.\n"
            "3) cumulative_kernels.h5, which contains the kernels used in the dose calculation.\n"
            "4) density.raw, which contains the raw density data of the phantom.\n"
            "5) omni_beam_lists.txt, which contains the specifications of the beam information"
            "All parameters should be specified as the arguments to this program")
        ("dataFolder", po::value<std::string>(), "The folder containing the necessary data (required)")
        ("resultFolder", po::value<std::string>(), "The folder to put result in (required)")
        ("deviceIdx", po::value<int>()->default_value(0), "The GPU idx for dose calculation. Here we only use 1 card.")
        ("unpack2Patient", po::value<bool>()->default_value(true), "The flag to unpack the BEV dose to the patient volume")
        ("debugLog", po::value<bool>()->default_value(false), "whether to log the variables")
        ("dicomVolumeDimension", po::value<std::vector<int>>()->multitoken(), 
            "Dicom volume dimension, 3 digits (required)")
        ("voxelSize", po::value<std::vector<float>>()->multitoken(),
            "Dicom voxel size. [cm] (required)")
        ("doseBoundingBoxStartIndices", po::value<std::vector<int>>()->multitoken(),
            "Dose bounding box start indices (required)")
        ("doseBoundingBoxDimensions", po::value<std::vector<int>>()->multitoken(),
            "Dose bounding box dimensions (required)")
        ("REVConvolutionArrayDimensions", po::value<std::vector<int>>()->multitoken(),
            "Maximum REV convolution array dimensions (required)")
        ("convlat", po::value<float>()->default_value(0.25),
            "Convolution ray lateral spacing. [cm]")
        ("convstep", po::value<float>()->default_value(0.25),
            "Convolution step spacing . [cm]")
        ("kernelExtent", po::value<float>()->default_value(4.0),
            "Dose kernel radius truncate distance. [cm]")
        ("nphi", po::value<int>()->default_value(8),
            "number of phi angles in convolution")
        ("ntheta",  po::value<int>()->default_value(8),
            "number of theta values in convolution")
        ("nradii", po::value<int>()->default_value(24),
            "number of radii values in convolution")
        ("penumbra", po::value<float>()->default_value(1.0),
            "Beamlet transverse dose spread. [cm]")
        ("beamCount", po::value<int>()->default_value(1),
            "Number of beams to pick from omni_beam_lists.txt");
        
        dicomVolumeStartCoords = std::vector<float>({0., 0., -25.55});
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if(vm.count("help"))
        {
            std::cout << desc << std::endl;
            return 1;
        }

        int width = 60;
        std::cout << "Parameters:" << std::endl;
        for (const auto& pair : vm)
        {
            std::stringstream second;
            const auto& value  = pair.second.value();
            if (auto ptr = boost::any_cast<int>(&value))
                second << *ptr;
            else if (auto ptr = boost::any_cast<float>(&value))
                second << *ptr;
            else if (auto ptr = boost::any_cast<std::vector<float>>(&value))
                second << *ptr;
            else if (auto ptr = boost::any_cast<std::vector<int>>(&value))
                second << *ptr;
            else if (auto ptr = boost::any_cast<std::string>(&value))
                second << *ptr;
            else
                second << "(unknown type)";
            
            std::string second_string = second.str();
            int remaining = width - pair.first.size() - second_string.size();
            remaining = std::max(5, remaining);

            std::stringstream output;
            output << pair.first << std::string(remaining, '.') << second_string;
            std::cout << output.str() << std::endl;
        }
        std::string first{"dicomVolumeStartCoords"};
        std::stringstream second;
        second << dicomVolumeStartCoords;
        std::string second_string = second.str();
        int remaining = width - first.size() - second_string.size();
        std::cout << first << std::string(remaining, '.') << second_string << std::endl;

        
        return 0;
}

const std::vector<float>& dev::getarg(const std::string& key)
{
    if (key == std::string("dicomVolumeStartCoords"))
            return dicomVolumeStartCoords;
        
    if (vm.count(key))
        return vm[key].as<std::vector<float>>();
    std::cerr << "The key " << key << " doesn't exist in the argument list!" << std::endl;
    exit(1);
}