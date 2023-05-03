#include <iostream>
#include "argparse.h"
#include "G4UIExecutive.hh"
#include "G4RunManagerFactory.hh"
#include "G4VisExecutive.hh"
#include "ActionInitialization.h"

#include "DetectorConstruction.h"
#include "QGSP_BERT.hh"

int main(int argc, char** argv)
{
    if (wk::argsInit(argc, argv)) return 0;
    
    G4UIExecutive* ui = nullptr;
    if (wk::getArg<bool>("gui"))
        ui = new G4UIExecutive(1, argv);
    
    auto* runManager = G4RunManagerFactory::CreateRunManager();

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    auto* detector = new wk::DetectorConstruction();
    runManager->SetUserInitialization(detector);

    G4VModularPhysicsList* physicsList = new QGSP_BERT;
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new wk::ActionInitialization);

    runManager->Initialize();
    
    if (ui)
    {
        ui->SessionStart();
        delete ui;
    }
    else
    {
        runManager->BeamOn(1);
    }
}
