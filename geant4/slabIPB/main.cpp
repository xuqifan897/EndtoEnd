#include "G4RunManagerFactory.hh"
#include "QGSP_BERT.hh"

#include "argparse.h"
#include "PhantomDef.h"
#include "DetectorConstruction.h"
#include "ActionInitialization.h"

int main(int argc, char** argv)
{
    if(si::argsInit(argc, argv))
        return 0;
    
    // geometry initialization
    si::GD = new si::GeomDef();
    
    auto* runManager = G4RunManagerFactory::CreateRunManager();

    // set random number seeds using time
    G4Random::setTheSeed(std::time(nullptr));

    auto* detector = new si::DetectorConstruction();
    runManager->SetUserInitialization(detector);

    G4VModularPhysicsList* physicsList = new QGSP_BERT;
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new si::ActionInitialization);

    runManager->Initialize();

    G4int nParticles = si::getArg<G4int>("nParticles");
    runManager->BeamOn(nParticles);

    delete runManager;
}