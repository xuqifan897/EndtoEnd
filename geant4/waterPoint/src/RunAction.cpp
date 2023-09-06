#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

#include "RunAction.h"
#include "Run.h"
#include "argparse.h"

G4Run* wp::RunAction::GenerateRun()
{
    return new Run();
}

void wp::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    if (this->isMaster)
        G4cout << "### Run: " << aRun->GetRunID() << "starts." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void wp::RunAction::EndOfRunAction(const G4Run* aRun)
{
    if (this->isMaster)
    {
        const Run* masterRun = static_cast<const Run*>(aRun);

        std::string resultFolder = (*vm)["resultFolder"].as<std::string>();
        fs::path ResFd(resultFolder);
        if (! fs::exists(ResFd))
            fs::create_directory(ResFd);

        int PhantomDimXY = (*vm)["PhantomDimXY"].as<int>();
        int PhantomDimZ = (*vm)["PhantomDimZ"].as<int>();
        int PhantomSZ = (*vm)["PhantomSZ"].as<int>();
        float Resolution = (*vm)["resolution"].as<float>() * cm;

        std::stringstream metadataSS;
        metadataSS << "Number of events: " << masterRun->GetNumberOfEvent() << std::endl;
        metadataSS << "Data type: double " << std::endl;
        metadataSS << "Dimension: (" << PhantomDimXY << ", " << PhantomDimXY << ", " 
            << PhantomDimZ << ") " << std::endl;
        metadataSS << "Source displacement: " << PhantomSZ << std::endl;
        metadataSS << "Resolution: " << Resolution / cm << " cm" << std::endl;
        metadataSS << "data order: z, y, x" << std::endl;

        std::string metadata = metadataSS.str();
        fs::path filePath = ResFd / fs::path("metadata.txt");
        std::ofstream file = std::ofstream(filePath.string());
        if (file.is_open())
        {
            file << metadata;
            file.close();
            std::cout << "Metadata written successfully!" << std::endl;
        }
        else
            std::cerr << "Cannot open file: " << filePath.string() << std::endl;

        std::vector<double> outputData(PhantomDimXY * PhantomDimXY * PhantomDimZ);
        auto HitsMaps = masterRun->getHitsMaps();
        int offset = 0;
        for (int i=0; i<HitsMaps.size(); i++)
        {
            const auto & map = *std::get<2>(HitsMaps[i])->GetMap();
            for (auto & it : map)
            {
                outputData[offset + it.first] = *(it.second);
            }
            offset += PhantomDimXY;
        }
        filePath = ResFd / fs::path("SD.bin");
        file = std::ofstream(filePath.string());
        if (file.is_open())
        {
            file.write((char*)(outputData.data()), outputData.size() * sizeof(double));
            file.close();
            std::cout << "Sensitive detector data written successfully!" << std::endl;
        }
        else
            std::cerr << "Unable to open file: " << filePath.string();
    }
}