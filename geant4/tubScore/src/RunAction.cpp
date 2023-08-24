#include "G4Run.hh"
#include "G4RunManager.hh"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

#include "RunAction.h"
#include "Run.h"
#include "argparse.h"
#include "PhantomDef.h"

G4Run* ts::RunAction::GenerateRun()
{
    return new Run();
}

void ts::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    if (this->isMaster)
        G4cout << "### Run: " << aRun->GetRunID() << "starts." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void ts::RunAction::EndOfRunAction(const G4Run* aRun)
{
    if (this->isMaster)
    {
        const Run* masterRun = static_cast<const Run*>(aRun);

        std::string resultFolder = (*ts::vm)["resultFolder"].as<std::string>();
        fs::path ResFd(resultFolder);
        if (! fs::exists(ResFd))
            fs::create_directory(ResFd);

        const auto & HitsMaps = masterRun->getHitsMaps();

        std::stringstream metadataSS;
        metadataSS << "Number of events: " << masterRun->GetNumberOfEvent() 
            << std::endl;
        metadataSS << "HitsMaps: " << std::endl;
        for (auto it=HitsMaps.begin(); it!=HitsMaps.end(); it++)
        {
            metadataSS << std::get<0>(*it) << "   dimension: "
                << std::get<2>(*it)->GetSize() << std::endl;
        }
        metadataSS << "Data type: double" << std::endl;
        std::string metadata = metadataSS.str();

        fs::path metaFile = ResFd / fs::path("metadata.txt");
        std::ofstream file(metaFile.string());
        if (file.is_open())
        {
            file << metadata;
            file.close();
            std::cout << "Metadata written successfully." << std::endl;
        }
        else
            std::cout << "Unable to open file: " << metaFile << std::endl;
        
        // prepare output data
        std::vector<double> output(GD->dimR * GD->dimZ);
        for (int i=0; i<HitsMaps.size(); i++)
        {
            const auto & hitsmap = *std::get<2>(HitsMaps[i]);
            int offset = i * GD->dimR;
            for (auto it : *(hitsmap.GetMap()))
                output[offset + it.first] = *(it.second);
        }

        std::string name = std::string("SD.bin");
        fs::path filePath = ResFd / fs::path(name);
        file = std::ofstream(filePath.string());
        if (file.is_open())
        {
            file.write((char*)(output.data()),
                output.size()*sizeof(double));
            file.close();
        }
        else
            std::cout << "Unable to open file: " << filePath.string() << std::endl;
    }
}