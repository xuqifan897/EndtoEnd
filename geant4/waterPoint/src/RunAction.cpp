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
        int PhantomBottom = (*vm)["PhantomBottom"].as<int>();
        float Resolution = (*vm)["resolution"].as<float>() * cm;

        std::stringstream metadataSS;
        metadataSS << "Number of events: " << masterRun->GetNumberOfEvent() << std::endl;
        metadataSS << "Data type: double " << std::endl;
        metadataSS << "Dimension: (" << PhantomDimXY << ", " << PhantomDimXY << ", " 
            << PhantomDimZ << ") " << std::endl;
        metadataSS << "Source displacement: " << PhantomSZ << std::endl;
        metadataSS << "Resolution: " << Resolution / cm << " cm" << std::endl;
        metadataSS << "data order: z, y, x" << std::endl;
    
        metadataSS << "Kernel Dimension: (" << PhantomDimXY << ", " 
            << PhantomDimXY << ", " << PhantomSZ + 
            PhantomBottom << ")" << std::endl;

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

        const auto & kernel = masterRun->GetKernel();
        filePath = ResFd / fs::path("SD.bin");
        file = std::ofstream(filePath.string());
        if (file.is_open())
        {
            for (int i=0; i<PhantomSZ + PhantomBottom; i++)
                for (int j=0; j<PhantomDimXY; j++)
                    file.write((char*)(kernel[i][j].data()), 
                        PhantomDimXY*sizeof(double));
            file.close();
            std::cout << "Kernel data written successfully!" << std::endl;
        }
        else
            std::cerr << "Unable to open file: " << filePath.string();
    }
}