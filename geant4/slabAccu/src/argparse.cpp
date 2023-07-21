#include <boost/program_options.hpp>
#include "globals.hh"

#include "argparse.h"

po::variables_map* sa::vm = nullptr;

int sa::argsInit(int argc, char** argv)
{
    po::options_description desc(
        "An application to realize accurate dose calculation in slab phantom");
    desc.add_options()
        ("help", "produce help messages. This program is to calculate the " 
            "infinitesimal pencil beam (IPB) dose distribution, in a slab "
            "phantom. The user can specify the geometry parameters in \""
            "./src/PhantomDef.cpp\". Resolution and sizes , "
            "are in half size, in accordance with the settings of Geant4")
        ("gui", po::value<bool>()->default_value(false),
            "Whether to use graphical user interface.")
        ("Energy", po::value<float>()->default_value(6.0),
            "the value of primary photon energy [MeV].")
        ("nParticles", po::value<int>()->default_value(1),
            "The number of particles to simulate.")
        
        // The following block is for logging
        ("resultFolder", po::value<std::string>(),
            "The name of the folder to which we write results");
    
    vm = new po::variables_map();
    po::store(po::parse_command_line(argc, argv, desc), *vm);
    po::notify(*vm);

    if ((*vm).count("help"))
    {
        G4cout << desc << G4endl;
        return 1;
    }

    G4cout << "Parameters:" << G4endl;
    for (const auto& pair : *vm)
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