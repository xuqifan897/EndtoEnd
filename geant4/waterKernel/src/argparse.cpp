#include "argparse.h"
#include <iostream>
#include <boost/program_options.hpp>

#include "G4ios.hh"

namespace po = boost::program_options;
po::variables_map wk::vm = 0;

int wk::argsInit(int argc, char** argv)
{
    po::options_description desc("The parameters for water kernel calculation");
    desc.add_options()
        ("help", "produce help messages")
        ("gui", po::value<bool>()->default_value(false), 
            "to use graphical user interface or not")
        ("energy", po::value<float>()->default_value(6.0), 
            "the value of primary photon energy [MeV]")
        ("sizeX", po::value<float>()->default_value(10.0), 
            "the phantom is a water box. This is the "
            "half length of the box along X axis [cm]")
        ("sizeY", po::value<float>()->default_value(10.0),
            "the phantom is a water box. This is the "
            "half length of the box along Y axis [cm]")
        ("sizeZ", po::value<float>()->default_value(10.0),
            "the phantom is a water box. This is the "
            "half length of the box along Z axis [cm]")
        ("posZ", po::value<float>()->default_value(-10.0),
            "the z component of the particle position. "
            "The X and Y position are 0. The particle "
            "shoots in the +z direction. [cm]")
        ("nParticles", po::value<int>()->default_value(1),
            "the number of particles to shoot. [int]")
        ("kernelType", po::value<std::string>()->default_value("IPB"), 
            "To select the type of dose kernel to calculate. "
            "One can choose either \"IPB\" or \"point\". The "
            "first refers to infinitesimal pencil beam kernel, "
            "the second refers to point kernel.\n"
            "In case of IPB, the scoring mesh is set to the "
            "same size as the world geometry, while in case "
            "of point kernel, the size should be specified by "
            "the user.")
        ("maxStep", po::value<float>()->default_value(-1.0), 
            "The maximum step size used in simulation, to allow "
            "for higher simulation accuracy. Does not apply if "
            "negative. [cm]")
        ("recordEventLog", po::value<bool>()->default_value(true),
            "To log the coordinates of individual hits.")
        ("kernelSizeX", po::value<float>()->default_value(10.0),
            "Only used in point kernel calculation. "
            "The half length of the scoring mesh in X dimension [cm]")
        ("kernelSizeY", po::value<float>()->default_value(10.0),
            "Only used in point kernel calculation. "
            "The half length of the scoring mesh in Y dimension [cm]")
        ("kernelSizeZ", po::value<float>()->default_value(10.0),
            "Only used in point kernel calculation. "
            "The half length of the scoring mesh in Z dimension [cm]")
        ("kernelPosZ", po::value<float>()->default_value(-5.0),
            "Only used in point kernel calculation. "
            "The z component of the interaction point w.r.t. the kernel. "
            "one should reserve some margin for back-scattering")
        ("kernelResX", po::value<float>()->default_value(0.05),
            "The resolution along the X dimension of either "
            "the IPB kernel or the point kernel, number of "
            "voxels would be calculated accordingly [cm]")
        ("kernelResY", po::value<float>()->default_value(0.05),
            "The resolution along the Y dimension of either "
            "the IPB kernel or the point kernel, number of "
            "voxels would be calculated accordingly [cm]")
        ("kernelResZ", po::value<float>()->default_value(0.05),
            "The resolution along the z dimension of either "
            "the IPB kernel or the point kernel, number of "
            "voxels would be calculated accordingly [cm]")
        ("kernelDimOdd", po::value<bool>()->default_value(true),
            "Whether the X and Y dimensions of the kernel to be odd")
        ("resultFolder", po::value<std::string>(),
            "The name of the folder to which we write results.")
        ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    // print options with values
    G4cout << "Parameters:" << G4endl;
    for (const auto& pair : vm)
    {
        G4cout << pair.first << " = ";
        const auto& value = pair.second.value();
        if (auto ptr = boost::any_cast<bool>(&value))
            G4cout << *ptr << G4endl;
        else if (auto ptr = boost::any_cast<int>(&value))
            G4cout << *ptr << G4endl;
        else if (auto ptr = boost::any_cast<float>(&value))
            G4cout << *ptr << G4endl;
        else if (auto ptr = boost::any_cast<std::string>(&value))
            G4cout << *ptr << G4endl;
        else
            G4cout << "(unknown type)" << G4endl;
    }

    return 0;
}