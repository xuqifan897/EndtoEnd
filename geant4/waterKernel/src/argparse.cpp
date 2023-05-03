#include "argparse.h"
#include <iostream>
#include <boost/program_options.hpp>

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
        ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    return 0;
}