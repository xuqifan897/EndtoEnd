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

#include "DetectorConstructionSI.h"
#include "PhantomDef.h"
#include "SensitiveDetector.h"
#include "config.h"
#include "argparse.h"

G4VPhysicalVolume* si::DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

#if PHANTOM == 4
    G4Element* elH = new G4Element(
        G4String("Hydrogen"), G4String("H"), 1., 1.01*g/mole);
    G4Element* elO = new G4Element(
        G4String("Ocygen"), G4String("O"), 8., 16.00*g/mole);

    G4double density = 1.000*g/cm3 * SPARSE;
    G4Material* water = new G4Material(G4String("water"), density, 2);
    water->AddElement(elH, 2);
    water->AddElement(elO, 1);

    density = 0.92*g/cm3 * SPARSE;
    G4Material* adipose = new G4Material(G4String("adipose"), density, 2);
    adipose->AddElement(elH, 2);
    adipose->AddElement(elO, 1);

    density = 1.04*g/cm3 * SPARSE;
    G4Material* muscle = new G4Material(G4String("muscle"), density, 2);
    muscle->AddElement(elH, 2);
    muscle->AddElement(elO, 1);

    density = 1.85*g/cm3 * SPARSE;
    G4Material* bone = new G4Material(G4String("bone"), density, 2);
    bone->AddElement(elH, 2);
    bone->AddElement(elO, 1);

    density = 0.25*g/cm3 * SPARSE;
    G4Material* lung = new G4Material(G4String("lung"), density, 2);
    lung->AddElement(elH, 2);
    lung->AddElement(elO, 1);

#else
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
    G4Material* adipose = nist->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");
    G4Material* muscle = nist->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP");
    G4Material* bone = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    G4Material* lung = nist->FindOrBuildMaterial("G4_LUNG_ICRP");
#endif
    std::map<G4String, G4Material*> LUTmat = {
        {"air", air},
        {"water", water},
        {"adipose", adipose},
        {"muscle", muscle},
        {"bone", bone},
        {"lung", lung}
    };

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

    G4bool checkOverlaps = true;

    auto solidWorld = new G4Box("World", GD->sizeX, GD->sizeY, GD->sizeZ);
    this->logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    auto physWorld = new G4PVPlacement(nullptr,
        G4ThreeVector(),
        logicWorld,
        "World",
        nullptr,
        false,
        0,
        checkOverlaps);
    
    // prepare user limits parameters
    G4double maxStep = si::getArg<float>("maxStep") * cm;

    this->logicals = std::vector<std::pair<G4String, G4LogicalVolume*>>();
    // prepare child slabs
    for (int i=0; i<GD->layers.size(); i++)
    {
        G4String layerName = G4String("Shape") + std::to_string(i+1);
        auto& layer = GD->layers[i];
        G4String& matName = std::get<0>(layer);
        G4double thickness = std::get<1>(layer);
        G4double offset = std::get<2>(layer);
        G4Material* material = LUTmat[matName];

        auto* solid = new G4Box(layerName, GD->sizeX, GD->sizeY, thickness);
        auto* logic = new G4LogicalVolume(solid, material, layerName);
        new G4PVPlacement(nullptr,
            G4ThreeVector(0., 0., offset),
            logic,
            layerName,
            logicWorld,
            false,
            0,
            checkOverlaps);
        this->logicals.push_back(std::make_pair(layerName, logic));

        //  set UserLimits
        if (maxStep > 0)
        {
            G4UserLimits* limit = new G4UserLimits(maxStep);
            logic->SetUserLimits(limit);
        }
    }
    return physWorld;
}

void si::DetectorConstruction::ConstructSDandField()
{

#if SENSDET == 0
    G4int index = 0;
    for (auto it=this->logicals.begin(); it!=this->logicals.end(); it++)
    {
        G4String& name = std::get<0>(*it);
        G4LogicalVolume* logical = std::get<1>(*it);
        si::SensitiveDetector* aSD = new SensitiveDetector(name, index);
        index ++;
        G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
        SetSensitiveDetector(logical, aSD);
    }
#elif SENSDET == 1
    G4String name("world");
    si::SensitiveDetector* aSD = new SensitiveDetector(name, 0);
    G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
    SetSensitiveDetector(this->logicWorld, aSD);

#elif SENSDET == 2
    G4String name("shared");
    si::SensitiveDetector* aSD = new SensitiveDetector(name, 0);
    G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
    for (auto it=this->logicals.begin(); it!=this->logicals.end(); it++)
    {
        G4LogicalVolume* logical = std::get<1>(*it);
        SetSensitiveDetector(logical, aSD);
    }
#endif
}