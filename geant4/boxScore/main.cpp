#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParallelWorldPhysics.hh"

#include "argparse.h"
#include "PhantomDef.h"
#include "DetectorConstruction.h"
#include "PhysicsList.h"
#include "ActionInitialization.h"
#include "Run.h"

int runGlobal(int argc, char** argv);

int main(int argc, char** argv)
{
    return runGlobal(argc, argv);
}

int runGlobal(int argc, char** argv)
{
    if (bs::argsInit(argc, argv))
        return 0;
    
    bs::GD = new bs::GeomDef();
    bs::GD->display();

    fs::path folder((*bs::vm)["resultFolder"].as<std::string>());
    if (! fs::exists(folder))
        fs::create_directories(folder);

    // prepare for the gloabl and local score matrices
    float voxelSize = (*bs::vm)["voxelSize"].as<float>() * cm;  // half voxel size
    float sizeZ = 0.;  // half size
    for (int i=0; i<bs::GD->layers.size(); i++)
        sizeZ += std::get<1>(bs::GD->layers[i]);
    int dimZ = std::round(sizeZ / voxelSize);
    std::cout << "Z dimension: " << dimZ;
    int dimXY = (*bs::vm)["dimXY"].as<int>();
    int SegZ = (*bs::vm)["SegZ"].as<int>();
    if (dimZ % SegZ != 0)
    {
        std::cerr << "dimZ is not a multiple of SegZ, error!" << std::endl;
        return 1;
    }
    std::vector<double> localScore(dimXY * dimXY * SegZ);

    G4Random::setTheSeed(std::time(nullptr));
    float offset = 0.;
    float thickness = SegZ * voxelSize;

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);

    auto detector = new bs::DetectorConstruction(offset, thickness);
    auto physicsList = new bs::PhysicsList();
    auto action = new bs::ActionInitialization(&localScore);

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(action);

    // write metadata
    std::ofstream metadataFile(folder / fs::path("metadata.txt"));
    metadataFile << "Block dimension (z, y, x): (" << SegZ << ", " << dimXY << ", " << dimXY << ")" << std::endl;
    metadataFile << "Number of blocks: " << dimZ / SegZ << std::endl;
    metadataFile << "Voxel size [cm] (half): " << voxelSize / cm << std::endl;
    metadataFile << "Data type: double" << std::endl;
    metadataFile.close();

    int nParticles = (*bs::vm)["nParticles"].as<int>();
    int iterations = dimZ / SegZ;
    for (int i=0; i<iterations; i++)
    {
        offset = i * SegZ * voxelSize;
        detector->getOffset() = offset;

        // reset score
        std::fill(localScore.begin(), localScore.end(), 0.);
        // reset count
        bs::Run::eventCounts.store(0);

        runManager->Initialize();
        runManager->BeamOn(nParticles);

        fs::path dataName = folder / fs::path("SD" + std::to_string(i+1) + ".bin");
        std::ofstream dataFile(dataName);
        if (dataFile.is_open())
        {
            dataFile.write((char*)(localScore.data()), SegZ*dimXY*dimXY*sizeof(double));
            dataFile.close();
            G4cout << "data " << i+1 << " written successfully!" << G4endl;
        }
        else
            G4cerr << "data " << i+1 << " writing unsuccessful." << G4endl;
    }

    delete runManager;
    return 0;
}