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

#include "DetectorConstruction.h"
#include "PhantomDef.h"
#include "SensitiveDetector.h"

si::DetectorConstruction::DetectorConstruction()
{}

si::DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* si::DetectorConstruction::Construct()
{
    // firstly, prepare the materials
    G4NistManager* NISTman = G4NistManager::Instance();
    G4Material* air = NISTman->FindOrBuildMaterial("G4_AIR");
    G4Material* water = NISTman->FindOrBuildMaterial("G4_WATER");
    G4Material* adipose = NISTman->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");
    G4Material* muscle = NISTman->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP");
    G4Material* bone = NISTman->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    G4Material* lung = NISTman->FindOrBuildMaterial("G4_LUNG_ICRP");
    std::map<G4String, G4Material*> LUTmat = {
        {"air", air},
        {"water", water},
        {"adipose", adipose},
        {"muscle", muscle},
        {"bone", bone},
        {"lung", lung}
    };

    // print materials
    for (auto it=LUTmat.begin(); it!=LUTmat.end(); it++)
        G4cout << it->second << G4endl;

    // get geometry properties
    if (GD == nullptr)
    {
        G4cout << "Material geometry is not initialized, "
            "please initialize it by calling \"si::GD = new "
            "si::GeomDef();\" before the detector construction" << G4endl;
        exit(1);
    }

    // construct world geometry
    G4Box* worldSolid = new G4Box(
        "world_solid", GD->sizeX, GD->sizeY, GD->sizeZ);
    G4LogicalVolume* worldLogical = new G4LogicalVolume(
        worldSolid, air, "world_logical", 0, 0, 0);
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(
        0,
        G4ThreeVector(),
        worldLogical,
        "world_physical",
        0,
        false,
        0,
        true);
    
    // construct child geometries
    G4int layerID = 1;
    for (auto it=GD->layers.begin(); it!=GD->layers.end(); it++)
    {
        G4String& material = std::get<0>(*it);
        G4double& thickness = std::get<1>(*it);
        G4double& offset = std::get<2>(*it);
        G4Material* Mat = LUTmat[material];

        G4String solidName = std::to_string(layerID) + material + "_solid";
        G4String logicalName = std::to_string(layerID) + material + "_logical";
        G4String physicalName = std::to_string(layerID) + material + "_physical";

        G4Box* solid = new G4Box(solidName, GD->sizeX, GD->sizeY, thickness);
        G4LogicalVolume* logical = new G4LogicalVolume(solid, Mat, logicalName, 0, 0, 0);
        G4PVPlacement(
            0,
            G4ThreeVector(0, 0, offset),
            logical,
            physicalName,
            worldLogical,
            false,
            0,
            true);
        
        this->logicals.push_back(std::make_pair(
            std::to_string(layerID)+material, logical));
        layerID ++;
    }
    return worldPhysical;
}


void si::DetectorConstruction::ConstructSDandField()
{
    for (auto it=this->logicals.begin(); it!=this->logicals.end(); it++)
    {
        G4String& name = std::get<0>(*it);
        G4LogicalVolume* logical = std::get<1>(*it);
        SensitiveDetector* aSD = new SensitiveDetector(name);
        G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
        SetSensitiveDetector(logical, aSD);
    }
}