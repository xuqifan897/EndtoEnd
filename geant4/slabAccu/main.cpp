#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "argparse.h"
#include "PhantomDef.h"
#include "DetectorConstruction.h"
#include "PhysicsList.h"
#include "ActionInitialization.h"

int main(int argc, char** argv)
{
    if (sa::argsInit(argc, argv))
        return 0;
    
    sa::GD = new sa::GeomDef;

    G4UIExecutive* ui = nullptr;
    if ((*sa::vm)["gui"].as<bool>())
        ui = new G4UIExecutive(argc, argv);
    
    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);
    
    // set random number seeds using time
    G4Random::setTheSeed(std::time(nullptr));

    runManager->SetUserInitialization(new sa::DetectorConstruction);

    runManager->SetUserInitialization(new sa::PhysicsList);

    runManager->SetUserInitialization(new sa::ActionInitialization);

    runManager->Initialize();

    int nParticles = (*sa::vm)["nParticles"].as<int>();
    runManager->BeamOn(nParticles);

    delete runManager;
}