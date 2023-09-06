#include <iostream>
#include "argparse.h"

po::variables_map* wp::vm=nullptr;

int wp::argsInit(int argc, char** argv)
{
    po::options_description desc(
        "An application to calculate the point dose kernel in water. "
    );
    desc.add_options()
        ("help", "An application to calculate the point dose kernel in water. "
        "Specifically, we shift the coordinates of the primary interaction "
        "to the origin, and collect the dose distribution")
        ("Energy", po::value<float>()->default_value(6.0),
            "The value of primary photon energy [MeV].")
        ("nParticles", po::value<int>()->default_value(1),
            "The number of particles to simulate")
        ("logFrequency", po::value<int>()->default_value(10000),
            "The frequency for logging events.")
        ("PrintTrajectory", po::value<bool>()->default_value(false),
            "Whether to print trajectory at the end of an event.")
        
        // phantom definition
        ("phantom", po::value<std::string>()->default_value("cube"),
            "can take the value of either cube or cylinder.")
        ("resolution", po::value<float>()->default_value(0.05), 
            "If using a Cartesian frame, this parameter means the "
            "half resolution along x, y, and z directions. If using a cylindrical "
            "frame, this parameter means the resolution along the radial "
            "direction and the half resolution along the z direction. [cm]")
        ("PhantomDimXY", po::value<int>()->default_value(200),
            "The phantom dimension along the x and y directions.")
        ("PhantomDimZ", po::value<int>()->default_value(400),
            "The phantom dimension along the z direction.")
        ("PhantomSZ", po::value<int>()->default_value(20),
            "The depth of the source. Assume the photon is along the +z "
            "direction, and is placed in the middle of x and y direction. It's "
            "a unitless argument.")
        ("TallyDimZ", po::value<int>()->default_value(200),
            "The tally dimension along the z direction. "
            "We set the tally X and Y dimensions to be the same as the phantom, "
            "and the source point is set to the same depth as phantom dim z.")
        
        // log information
        ("resultFolder", po::value<std::string>(), "The folder to which we write results");

    vm = new po::variables_map();
    po::store(po::parse_command_line(argc, argv, desc), *vm);
    po::notify(*vm);

    if ((*vm).count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    std::cout << "Parameters:" << std::endl;
    for (const auto& pair : *vm)
    {
        std::cout << pair.first << " = ";
        const auto& value = pair.second.value();
        if (auto ptr = boost::any_cast<bool>(&value))
            std::cout << *ptr << std::endl;
        else if (auto ptr = boost::any_cast<int>(&value))
            std::cout << *ptr << std::endl;
        else if (auto ptr = boost::any_cast<float>(&value))
            std::cout << *ptr << std::endl;
        else if (auto ptr = boost::any_cast<std::string>(&value))
            std::cout << *ptr << std::endl;
        else
            std::cout << "(unknown type)" << std::endl;
    }
    std::cout << std::endl;

    return 0;
}