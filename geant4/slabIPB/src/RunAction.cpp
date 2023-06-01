#include <fstream>
#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#include "argparse.h"
#include "RunAction.h"
#include "Run.h"

G4Run* si::RunAction::GenerateRun()
{
    return new si::Run();
}

void si::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    if (IsMaster())
        G4cout << "### Run: " << aRun->GetRunID() << " starts." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void si::RunAction::EndOfRunAction(const G4Run* aRun)
{
    if (IsMaster())
    {
        const si::Run* masterRun = static_cast<const si::Run*>(aRun);

        G4String resultFolder = getArg<std::string>("resultFolder");
        fs::path ResFd(resultFolder);
        if (! fs::exists(ResFd))
            fs::create_directory(ResFd);

        std::stringstream metadataSS;
        metadataSS << "Number of events: " << masterRun->GetNumberOfEvent() << G4endl;
        metadataSS << "Hits map dimension: (" << masterRun->getDimX() << ", "
            << masterRun->getDimY() << ", " << masterRun->getDimZ() << ")" << G4endl;
        metadataSS << "Hits map resolution: (" << masterRun->getFullResX()/cm << ", "
            << masterRun->getFullResY()/cm << ", " << masterRun->getFullResZ()/cm 
            << ")[cm]" << G4endl;
        metadataSS << "Data type: double" << G4endl;
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
        masterRun->writeHitsMap(arrayData.string());
    }
}