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
#include "ParallelWorld.h"
#include "PhysicsList.h"
#include "ActionInitialization.h"

int main(int argc, char** argv)
{
    if (bs::argsInit(argc, argv))
        return 0;
    
    G4UIExecutive* ui = nullptr;
    ui = new G4UIExecutive(argc, argv);
    
    bs::GD = new bs::GeomDef();
    bs::GD->display();

    // set random number seeds using time

    // prepare for the gloabl and local score matrices
    float voxelSize = (*bs::vm)["voxelSize"].as<float>() * cm;  // half voxel size
    float sizeZ = 0.;  // half size
    for (int i=0; i<bs::GD->layers.size(); i++)
        sizeZ += std::get<1>(bs::GD->layers[i]);
    int dimZ = std::round(sizeZ / voxelSize);
    std::cout << "Z dimension: " << dimZ;
    int dimXY = (*bs::vm)["dimXY"].as<int>();
    int SegZ = (*bs::vm)["SegZ"].as<int>();
    // std::vector<double> globalScore(dimXY * dimXY * dimZ);
    std::vector<double> localScore(dimXY * dimXY * SegZ);


    G4Random::setTheSeed(std::time(nullptr));
    G4String parallelWorldName = "ReadoutWorld";
    float offset = 0.;
    G4VUserDetectorConstruction* detector = new bs::DetectorConstruction();
    if ((*bs::vm)["scoring"].as<bool>())
        detector->RegisterParallelWorld(new bs::ParallelWorld(parallelWorldName, offset, SegZ));

    G4VModularPhysicsList* physicsList = new bs::PhysicsList();
    if ((*bs::vm)["scoring"].as<bool>())
        physicsList->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);

    runManager->SetUserInitialization(detector);

    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new bs::ActionInitialization(&localScore));

    runManager->Initialize();

    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    if(ui)
    {
        G4String command = "/control/execute ./vis.mac";
        G4UImanager::GetUIpointer()->ApplyCommand(command);
        ui->SessionStart();
        delete ui;
    }
    else
    {
        int nParticles = (*bs::vm)["nParticles"].as<int>();
        runManager->BeamOn(nParticles);
    }
    
    if (visManager)
        delete visManager;
    delete runManager;


    // data io
    fs::path folder((*bs::vm)["resultFolder"].as<std::string>());
    if (! fs::exists(folder))
        fs::create_directories(folder);
    fs::path file = folder / fs::path("SD.bin");
    std::ofstream File(file.string());
    if (File.is_open())
    {
        File.write((char*)(localScore.data()), localScore.size()*sizeof(double));
        File.close();
    }
    else
        std::cerr << "Unable to open file: " << file.string() << std::endl;

    std::stringstream log;
    log << "dimension (z, y, x) = (" << SegZ << ", " << 
        dimXY << ", " << dimXY << ")" << std::endl;
    log << "data type : double" << std::endl;
    log << "offset = " << offset / cm << "cm" << std::endl;
    file = folder / fs::path("metadata.txt");
    File = std::ofstream(file.string());
    if (File.is_open())
    {
        File << log.str();
        File.close();
    }
    else
        std::cerr << "Unable to open file: " << file.string() << " " << std::endl;
}