#include "ParallelWorld.h"
#include "argparse.h"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"

bs::ParallelWorld::ParallelWorld(G4String& parallelWorldName, float offset, int dimz)
: G4VUserParallelWorld(parallelWorldName), Offset(offset), DimZ(dimz)
{
    this->DimXY = (*bs::vm)["dimXY"].as<int>();
    this->voxelSize = (*bs::vm)["voxelSize"].as<float>() * cm;
}

void bs::ParallelWorld::Construct()
{
    G4Material* dummyMat = nullptr;

    G4VPhysicalVolume* WorldPhys = GetWorld();
    G4LogicalVolume* WorldLog = WorldPhys->GetLogicalVolume();
    float SizeXY = this->DimXY * this->voxelSize;
    bool checkOverlaps = false;
    this->logicals.reserve(this->DimZ * this->DimXY);

    for (int i=0; i<this->DimZ; i++)
    {
        // iterate over dim z
        float offsetZ = this->Offset + (2 * i + 1) * this->voxelSize;
        for (int j=0; j<this->DimXY; j++)
        {
            // iterate over dim y
            // Create a bar of detectors
            std::string barName = "bar_" + std::to_string(i+1) + 
                "_" + std::to_string(j+1);
            auto barS = new G4Box(barName, SizeXY, this->voxelSize, this->voxelSize);
            auto barLV = new G4LogicalVolume(
                barS,
                dummyMat,
                barName);
            float offsetY = - SizeXY + (2 * j + 1) * this->voxelSize;
            new G4PVPlacement(
                nullptr,  // no rotation
                G4ThreeVector(0, offsetY, offsetZ),  // offset
                barLV,  // its logical volume
                barName,  // its name
                WorldLog,  // its mother volume
                false,  //  no boolean operator
                0,  // copy number
                checkOverlaps);
            
            // subdivide detector elements through X dimension
            std::string barElementName = "barElement_" + std::to_string(i+1) + 
                std::to_string(j+1);
            auto elementS = new G4Box(barElementName, this->voxelSize, 
                this->voxelSize, this->voxelSize);
            auto elementLV = new G4LogicalVolume(elementS, dummyMat, barElementName);
            new G4PVReplica(
                barElementName,  // detector element name
                elementLV,  // its logical volume
                barLV,  // mother logical volume
                kXAxis,  // axis
                this->DimXY,  // n replicas
                2*this->voxelSize);  // width
            this->logicals.push_back(elementLV);
        }
    }
}

void bs::ParallelWorld::ConstructSD()
{
    auto SDMpointer = G4SDManager::GetSDMpointer();
    SDMpointer->SetVerboseLevel(1);
    
    for (int i=0; i<this->logicals.size(); i++)
    {
        std::string SDname = std::string("SD") + std::to_string(i+1);
        auto senseDet = new G4MultiFunctionalDetector(SDname);
        SDMpointer->AddNewDetector(senseDet);

        G4VPrimitiveScorer* primitive;
        primitive = new G4PSEnergyDeposit("Edep");
        senseDet->RegisterPrimitive(primitive);

        SetSensitiveDetector(this->logicals[i], senseDet);
    }
}