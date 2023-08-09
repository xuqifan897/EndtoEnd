#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "argparse.h"

int main(int argc, char** argv)
{
    if (ts::argsInit(argc, argv))
        return 0;
    
    
}