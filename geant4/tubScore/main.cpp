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
    if (ts::argsInit(argc, argv))
        return 0;
    
    ts::GD = new ts::GeomDef();

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);

    // set random number seeds using time
    G4Random::setTheSeed(std::time(nullptr));

    runManager->SetUserInitialization(new ts::DetectorConstruction);

    runManager->SetUserInitialization(new ts::PhysicsList);

    runManager->SetUserInitialization(new ts::ActionInitialization);

    runManager->Initialize();

    int nParticles = (*ts::vm)["nParticles"].as<int>();
    runManager->BeamOn(nParticles);

    delete runManager;
}