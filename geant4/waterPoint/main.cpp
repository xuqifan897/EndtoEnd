// This project aims to calculate the point dose kernel in 
// water phantom, to better formulate the dose distribution 
// algorithm in inhomogeneous phantoms

#include "G4RunManagerFactory.hh"

#include "argparse.h"
#include "DetectorConstruction.h"
#include "PhysicsList.h"
#include "ActionInitialization.h"

int main(int argc, char** argv)
{
    if (wp::argsInit(argc, argv))
        return 0;
    
    auto* runManager = G4RunManagerFactory::CreateRunManager();

    G4Random::setTheSeed(std::time(nullptr));

    runManager->SetUserInitialization(new wp::DetectorConstruction);

    runManager->SetUserInitialization(new wp::PhysicsList);

    runManager->SetUserInitialization(new wp::ActionInitialization);

    runManager->Initialize();

    int nParticles = (*wp::vm)["nParticles"].as<int>();
    runManager->BeamOn(nParticles);

    delete runManager;
}