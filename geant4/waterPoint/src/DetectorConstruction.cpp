#include "iostream"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4Threading.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4VUserDetectorConstruction.hh"

#include "DetectorConstruction.h"
#include "argparse.h"


G4VPhysicalVolume* wp::DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

    std::cout << water;
    std::string phantomType = (*vm)["phantom"].as<std::string>();
    int PhantomDimXY = (*vm)["PhantomDimXY"].as<int>();
    int PhantomDimZ = (*vm)["PhantomDimZ"].as<int>();
    float resolution = (*vm)["resolution"].as<float>() * cm;
    bool checkOverlaps = true;

    G4VPhysicalVolume* worldPV = nullptr;

    if (phantomType == std::string("cube"))
    {
        std::cout << "A cube phantom is used." << std::endl;
        auto worldSV = new G4Box("world", PhantomDimXY*resolution, 
            PhantomDimXY*resolution, PhantomDimZ*resolution);
        
        auto worldLV = new G4LogicalVolume(worldSV, water, "World");
        worldPV = new G4PVPlacement(
            nullptr,  // no rotation
            G4ThreeVector(),  // at (0, 0, 0)
            worldLV,  // its logical volume
            "World",  // its name
            nullptr,  // no mother volume
            false,  // no boolean operator
            0,  // copy number
            checkOverlaps);

        this->SenseDetList = std::vector<G4LogicalVolume*>(PhantomDimZ * PhantomDimXY);
        // order z, y, x
        for (int i=0; i<PhantomDimZ; i++)
        {
            float offsetZ = (2 * i + 1 - PhantomDimZ) * resolution;
            for (int j=0; j<PhantomDimXY; j++)
            {
                std::string name = std::string("slabX_") + std::to_string(i+1) + 
                    std::string("_") + std::to_string(j+1);
                auto solidX = new G4Box(name, PhantomDimXY*resolution, resolution, resolution);
                auto logicalX = new G4LogicalVolume(solidX, water, name);
                float offsetY = (2 * j + 1 - PhantomDimXY) * resolution;
                new G4PVPlacement(
                    nullptr,  // no rotation
                    G4ThreeVector(0., offsetY, offsetZ),  // offset
                    logicalX,  // its logical volume
                    name,  // its name
                    worldLV,  // mother volume
                    false,  // no boolean operator
                    0,  // copy number
                    checkOverlaps
                );
                // Repeat along the x axis
                name = std::string("slabXX_") + std::to_string(i+1) + 
                    std::string("_") + std::to_string(j+1);
                auto solidXX = new G4Box(name, resolution, resolution, resolution);
                auto logicalXX = new G4LogicalVolume(solidXX, water, name);
                new G4PVReplica(name, logicalXX, logicalX, kXAxis, PhantomDimXY, 2*resolution);

                int idx = i*PhantomDimXY + j;
                this->SenseDetList[idx] = logicalXX;
            }
        }
    }
    else if (phantomType == std::string("cylinder"))
    {
        std::cout << "A cylinder phantom is used." << std::endl;
    }
    else
    {
        std::cerr << "unrecognized phantom type: " << phantomType << std::endl;
        exit(EXIT_FAILURE);
    }

    return worldPV;
}

void wp::DetectorConstruction::ConstructSDandField()
{
    //  Sensitive Detector Manager.
    G4SDManager* pSDman = G4SDManager::GetSDMpointer();

    std::string phantomType = (*vm)["phantom"].as<std::string>();
    if (phantomType == "cube")
    {
        for (int i=0; i<this->SenseDetList.size(); i++)
        {
            G4String SDname = std::string("det") + std::to_string(i+1);
            G4MultiFunctionalDetector* det = 
                new G4MultiFunctionalDetector(SDname);
            pSDman->AddNewDetector(det);

            auto* primitive = new G4PSEnergyDeposit("Edep");
            det->RegisterPrimitive(primitive);

            SetSensitiveDetector(this->SenseDetList[i], det);
        }
    }
    else if (phantomType == "cylinder")
    {
        std::cerr << "Cylinder sensitive detector construction undefined." << std::endl;
        exit(EXIT_FAILURE);
    }
    return;
}