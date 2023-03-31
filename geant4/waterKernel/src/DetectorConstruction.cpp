#include "DetectorConstruction.h"
#include "argparse.h"

#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

wk::DetectorConstruction::DetectorConstruction()
    :   G4VUserDetectorConstruction()
{
    sizeX = getArg<float>("sizeX") * cm;
    sizeY = getArg<float>("sizeY") * cm;
    sizeZ = getArg<float>("sizeZ") * cm;
}

wk::DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* wk::DetectorConstruction::Construct()
{
    // firstly, prepare the material
    G4NistManager* NISTman = G4NistManager::Instance();
    G4Material* water = NISTman->FindOrBuildMaterial("G4_WATER");

    G4Box* experimentalHall_box = 
        new G4Box("expHall_b", sizeX, sizeY, sizeZ);
    experimentalHall_log = new G4LogicalVolume(
        experimentalHall_box, water, "expHall_L", 0, 0, 0);
    G4VPhysicalVolume* experimentalHall_phys = new G4PVPlacement(
        0, G4ThreeVector(), experimentalHall_log, "expHall_P", 0, false, 0);

    
    return experimentalHall_phys;
}