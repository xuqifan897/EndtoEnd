#include "G4RunManagerFactory.hh"
#include "QGSP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

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

    G4UIExecutive* ui = nullptr;
    G4VisManager* visManager = nullptr;
    G4UImanager* UImanager = nullptr;
    G4bool gui = si::getArg<bool>("gui");
    if (gui)
    {
        ui = new G4UIExecutive(1, argv);
        visManager = new G4VisExecutive;
        visManager->Initialize();
        UImanager = G4UImanager::GetUIpointer();
    }

    auto* runManager = G4RunManagerFactory::CreateRunManager();

    // set random number seeds using time
    G4Random::setTheSeed(std::time(nullptr));

    auto* detector = new si::DetectorConstruction();
    runManager->SetUserInitialization(detector);

    G4VModularPhysicsList* physicsList = new QGSP_BERT;
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new si::ActionInitialization);

    runManager->Initialize();

    if (gui)
    {
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }
    else
    {
        G4int nParticles = si::getArg<G4int>("nParticles");
        runManager->BeamOn(nParticles);
    }

    delete visManager;
    delete runManager;
}