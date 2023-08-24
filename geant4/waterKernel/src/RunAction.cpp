#include "RunAction.h"
#include "wkRun.h"
#include "config.h"

#include <fstream>
#include <string>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

wk::RunAction::RunAction()
    :G4UserRunAction()
{}

wk::RunAction::~RunAction()
{}

#if PHASE > 0
G4Run* wk::RunAction::GenerateRun()
{ return new wk::Run; }
#endif

void wk::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

void wk::RunAction::EndOfRunAction(const G4Run* aRun)
{
#if PHASE > 0
    const Run* bRun = static_cast<const Run* >(aRun);
    if (IsMaster())
    {
        G4cout << "Number of valid events (whose first interaction "
            "points are within the valid range): " << bRun->GetNumberOfEvent() << G4endl;
        
        G4String resultFolder = getArg<std::string>("resultFolder");
        fs::path ResFd(resultFolder);
        if (! fs::exists(ResFd))
            fs::create_directories(ResFd);
        
        std::stringstream metadataSS;
        metadataSS << "Number of valid events (whose first interaction "
            "points are within the valid range): " << bRun->GetNumberOfEvent() 
            << std::endl;
        metadataSS << "Hits map dimension: (" << bRun->getDimX() << ", " 
            << bRun->getDimY() << ", " << bRun->getDimZ() << ")" << std::endl;
        metadataSS << "Hits map resolution: (" << bRun->getFullResX()/cm << ", "
            << bRun->getFullResY()/cm << ", " << bRun->getFullResZ()/cm 
            << ")[cm]" << std::endl;
        metadataSS << "Kernel source offset: (" << bRun->getPosZ()/cm 
            << ")[cm]" << std::endl;
        metadataSS << "Data type: double";
        std::string metadata = metadataSS.str();

        fs::path metaFile = ResFd / fs::path("metadata.txt");
        std::ofstream file(metaFile.string());
        if (file.is_open())
        {
            file << metadata;
            file.close();
            G4cout << "Metadata written successfully." << G4endl;
        }
        else
            G4cout << "Unable to open file: " << metaFile << G4endl;
        
        fs::path arrayData = ResFd / fs::path("array.bin");
        bRun->writeHitsMap(arrayData.string());
    }
#endif
}