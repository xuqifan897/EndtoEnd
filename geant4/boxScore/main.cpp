#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;
#include <string>

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
#include "AlternaltiveDetectorConstruction.h"
#include "DummyActionInitialization.h"

// int runGlobal(int argc, char** argv);
int runLocal(int argc, char** argv);

int main(int argc, char** argv)
{
    return runLocal(argc, argv);
}


#if false
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
#endif

int runLocal(int argc, char** argv)
{
    if (bs::argsInit(argc, argv))
        return 0;
    bool gui = (*bs::vm)["gui"].as<bool>();

    G4UIExecutive* ui = nullptr;
    if (gui)
        ui = new G4UIExecutive(argc, argv);

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);
    
    G4VisManager* visManager = nullptr;
    if (gui)
    {
        visManager = new G4VisExecutive;
        visManager->Initialize();
    }
    
    bs::GD = new bs::GeomDef();
    bs::GD->display();

    fs::path folder((*bs::vm)["resultFolder"].as<std::string>());
    if (! fs::exists(folder))
        fs::create_directories(folder);

    // prepare for the gloabl and local score matrices
    float voxelSize = (*bs::vm)["voxelSize"].as<float>() * cm;  // half voxel size
    int iteration = (*bs::vm)["iteration"].as<int>();  // the iteration of this execution
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
    if (iteration >= dimZ / SegZ)
    {
        std::cerr << "iteration: " << iteration << ", dimZ / SegZ: " << dimZ / SegZ << std::endl;
        return 1;
    }
    std::vector<double> localScore(dimXY * dimXY * SegZ);

    // write metadata
    std::stringstream metadataName;
    metadataName << "metadata_" << std::setw(3) << std::setfill('0') << iteration+1 << ".txt";
    std::ofstream metadataFile(folder / fs::path(metadataName.str()));
    metadataFile << "Block dimension (z, y, x): (" << SegZ << ", " << dimXY << ", " << dimXY << ")" << std::endl;
    metadataFile << "Number of blocks: " << dimZ / SegZ << std::endl;
    metadataFile << "Voxel size [cm] (half): " << voxelSize / cm << std::endl;
    metadataFile << "Data type: double" << std::endl;
    metadataFile.close();

    G4Random::setTheSeed(std::time(nullptr));
    float thickness = SegZ * voxelSize;
    float offset = iteration * thickness;

    G4VUserDetectorConstruction* detector = nullptr;
    G4VUserActionInitialization* action = nullptr;
    if ((*bs::vm)["dummy"].as<bool>())
    {
        detector = new bs::AlternativeDetectorConstruction(offset, thickness);
        action = new bs::DummyActionInitialization;
    }
    else
    {
        detector = new bs::DetectorConstruction(offset, thickness);
        action = new bs::ActionInitialization(&localScore);
    }
    auto physicsList = new bs::PhysicsList();

    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(action);
    runManager->Initialize();

    if (gui)
    {
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        std::string command("/control/execute init_vis.mac");
        UImanager->ApplyCommand(command);
        ui->SessionStart();
        
        delete ui;
        delete visManager;
        delete runManager;
        return 1;
    }

    int nParticles = (*bs::vm)["nParticles"].as<int>();
    runManager->BeamOn(nParticles);

    // log data
    std::stringstream dataNameSS;
    dataNameSS << "SD" << std::setw(3) << std::setfill('0') << iteration+1 << ".bin";
    fs::path dataName = folder / fs::path(dataNameSS.str());
    std::ofstream dataFile(dataName);
    if (dataFile.is_open())
    {
        dataFile.write((char*)(localScore.data()), SegZ*dimXY*dimXY*sizeof(double));
        dataFile.close();
        G4cout << "data " << iteration+1 << " written successfully!" << G4endl;
    }
    else
        G4cerr << "data " << iteration+1 << " writing unsuccessful." << G4endl;

    delete runManager;
    return 0;
}