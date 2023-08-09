#include <boost/program_options.hpp>
#include "globals.hh"
#include <string>
#include <iostream>

#include "argparse.h"

po::variables_map* ts::vm = nullptr;

int ts::argsInit(int argc, char** argv)
{
    po::options_description desc(
        "An application to calculate both the longitudinal "
        "and transversal dose distribution");
    desc.add_options()
        ("help", "This program is to calculate both the longitudinal "
            "and transversal dose distribution of an infinitesimal "
            "pencil beam (IPB), in a slab phantom. The user can specify "
            "the geometry in ./src/PhantomDef.cpp. Resolution and sizes "
            "are in half, in accordance with the settings of Geant4")
        ("Energy", po::value<float>()->default_value(6.),
            "the value of primary photon energy [MeV].")
        ("nParticles", po::value<int>()->default_value(1),
            "The number of particles to simulate")
        
        // The following block is for logging
        ("resultFolder", po::value<std::string>(),
            "The folder to which we write results");
    
    vm = new po::variables_map();
    po::store(po::parse_command_line(argc, argv, desc), *vm);
    po::notify(*vm);

    if ((*vm).count("help"))
    {
        std::cout << desc << std::endl;
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