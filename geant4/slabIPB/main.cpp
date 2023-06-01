#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "argparse.h"
#include "DetectorConstructionB1.h"
#include "DetectorConstructionSI.h"
#include "ActionInitialization.h"
#include "PhantomDef.h"

int main(int argc, char** argv)
{
    if (si::argsInit(argc, argv))
        return 0;
    
    si::GD = new si::GeomDef;

    G4UIExecutive* ui = nullptr;
    if (si::getArg<bool>("gui"))
        ui = new G4UIExecutive(1, argv);

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);
    
    // set random number seeds using time
    G4Random::setTheSeed(std::time(nullptr));
    
    // runManager->SetUserInitialization(new si::DetectorConstructionB1());
    runManager->SetUserInitialization(new si::DetectorConstruction);

    G4VModularPhysicsList* physicsList = new QBBC;
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new si::ActionInitialization);
    runManager->Initialize();

    G4VisManager* visManager = nullptr;
    G4UImanager* UImanager = nullptr;
    if (si::getArg<bool>("gui"))
    {
        visManager = new G4VisExecutive;
        visManager->Initialize();

        UImanager = G4UImanager::GetUIpointer();

        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
        delete visManager;
    }
    else
    {
        G4int nParticles = si::getArg<G4int>("nParticles");
        runManager->BeamOn(nParticles);
    }
    
    delete runManager;
}