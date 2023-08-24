#include "DetectorConstruction.h"
#include "argparse.h"
#include "TrackerSD.h"

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
        new G4Box("expHall_b", this->sizeX, this->sizeY, this->sizeZ);
    this->experimentalHall_log = new G4LogicalVolume(
        experimentalHall_box, water, "expHall_L", 0, 0, 0);
    G4VPhysicalVolume* experimentalHall_phys = new G4PVPlacement(
        0, G4ThreeVector(), this->experimentalHall_log, "expHall_P", 0, false, 0);

    // set some attributes
    G4VisAttributes* experimentalHallVisAtt = 
        new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    experimentalHallVisAtt->SetForceWireframe(true);
    this->experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);

    // try to set maxStep for higher-accuracy simulation
    G4double maxStep = getArg<float>("maxStep") * cm;
    if (maxStep < 0)
        G4cout << "Do not set maxStep." << G4endl;
    else
    {
        G4UserLimits* userLimits = new G4UserLimits();
        userLimits->SetMaxAllowedStep(maxStep);
        this->experimentalHall_log->SetUserLimits(userLimits);
    }
    
    return experimentalHall_phys;
}

void wk::DetectorConstruction::ConstructSDandField()
{
    // define sensitive detector
    G4String trackerSDname = "mydet";
    TrackerSD* trackerSD = new TrackerSD(trackerSDname);
    G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD);
    SetSensitiveDetector(this->experimentalHall_log, trackerSD);
}