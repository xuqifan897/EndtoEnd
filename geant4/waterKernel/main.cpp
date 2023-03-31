#include <iostream>
#include "argparse.h"
#include "G4UIExecutive.hh"
#include "G4RunManagerFactory.hh"
#include "G4VisExecutive.hh"

int main(int argc, char** argv)
{
    if (wk::argsInit(argc, argv)) return 0;
    
    G4UIExecutive* ui = nullptr;
    if (wk::getArg<bool>("gui"))
        ui = new G4UIExecutive(1, argv);
    
    auto* runManager = G4RunManagerFactory::CreateRunManager();

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    
}
