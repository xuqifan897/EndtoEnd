#include <boost/program_options.hpp>
#include "globals.hh"
namespace po = boost::program_options;

#include "argparse.h"

po::variables_map si::vm = 0;

int si::argsInit(int argc, char** argv)
{
    po::options_description desc("The parameters for IPB dose calculation");
    desc.add_options()
        ("help", "produce help messages. This program is to calculate the "
            "infinitesimal pencil beam (IPB) dose distribution, for either water "
            "or slab phantom. To switch the type of phantom, the user should "
            "uncomment either line 10 or 11 of CMakeLists.txt. To change the "
            "geometry parameters, the user should change the parameters in "
            "\"./src/PhantomDef.cpp.\" All dimensions, including resolution, "
            "are in half size, in accordance with the settings of Geant4.")
        ("gui", po::value<bool>()->default_value(false),
            "Whether to use graphical user interface ")
        ("Energy", po::value<float>()->default_value(6.0),
            "the value of primary photon energy [MeV]")
        ("nParticles", po::value<G4int>()->default_value(1),
            "The number of particles to simulate.")
        ("sourceZ", po::value<float>()->default_value(-15.0),
            "The Z position of the source. [cm]. This entry was only used in "
            "the test phase. In real applications, it is determined by the "
            "geometry setting")
        
        // The following block is for UserLimits
        ("maxStep", po::value<float>()->default_value(-1),
            "The maximum step size allowed in the simulation. "
            "The lower the step size, the higher the simulation "
            "accuracy. When it's negative, it's not enforced. [cm]")
        ("maxTrack", po::value<float>()->default_value(-1),
            "The maximum allowed track length. When negative, not enforced. [cm]")
        ("maxTime", po::value<float>()->default_value(-1),
            "Maximum allowed time per step. When negative, not enforced [ns]")
        ("minEkine", po::value<float>()->default_value(0),
            "The minimum energy cut for stacking secondaries. [MeV]")
        ("minRange", po::value<float>()->default_value(0),
            "The minimum range cut for stacking secondaries. "
            "The same function as above, but translates to different "
            "energies in different materials. [cm]")

        // The following block is for logging
        ("resultFolder", po::value<std::string>(),
            "The name of the folder to which we write results.")
        ("recordEventLog", po::value<bool>()->default_value(false),
            "To log the coordinates of individual hits.");
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        G4cout << desc << G4endl;
        return 1;
    }

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