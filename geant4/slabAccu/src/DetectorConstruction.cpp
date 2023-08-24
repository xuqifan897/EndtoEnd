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
#include "PhantomDef.h"
#include "config.h"

G4VPhysicalVolume* sa::DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

    G4Element* elH = new G4Element(
        std::string("Hydrogen"), std::string("H"), 1., 1.01*g/mole);
    G4Element* elO = new G4Element(
        std::string("Ocygen"), std::string("O"), 8., 16.00*g/mole);

    G4double density = 1.000*g/cm3;
    G4Material* water = new G4Material(std::string("water"), density, 2);
    water->AddElement(elH, 2);
    water->AddElement(elO, 1);

    density = 0.92*g/cm3;
    G4Material* adipose = new G4Material(std::string("adipose"), density, 2);
    adipose->AddElement(elH, 2);
    adipose->AddElement(elO, 1);

    density = 1.04*g/cm3;
    G4Material* muscle = new G4Material(std::string("muscle"), density, 2);
    muscle->AddElement(elH, 2);
    muscle->AddElement(elO, 1);

    density = 1.85*g/cm3;
    G4Material* bone = new G4Material(std::string("bone"), density, 2);
    bone->AddElement(elH, 2);
    bone->AddElement(elO, 1);

    density = 0.25*g/cm3;
    G4Material* lung = new G4Material(std::string("lung"), density, 2);
    lung->AddElement(elH, 2);
    lung->AddElement(elO, 1);

    std::map<std::string, G4Material*> LUTmat = {
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

    bool checkOverlaps = true;

    // firstly, construct the world volume
    auto worldS = new G4Box("World",  // name
        sa::GD->sizeX, sa::GD->sizeY, sa::GD->sizeZ);  // size

    auto worldLV = new G4LogicalVolume(
        worldS,  // solid
        air,  // material
        "World");  // name

    auto worldPV = new G4PVPlacement(nullptr,  // no rotation
        G4ThreeVector(),  // at (0, 0, 0)
        worldLV,  // its logical volume
        "World",  // its name
        nullptr,  // no mother volume
        false,    // no boolean operator
        0,        // copy number
        checkOverlaps);  // checking overlaps

#if PHANTOM == 1
    // for this example, we are going to test the accuracy of the native 
    // G4MultiFunctionalDetector class. In the layers.

    for (int i=0; i<sa::GD->layers.size(); i++)
    {
        std::string name = std::string("layer") + std::to_string(i+1);
        const auto & param = GD->layers[i];
        const auto & thickness = std::get<1>(param);
        const auto & displacement = std::get<2>(param);

        auto layerS = new G4Box(name,
            sa::GD->sizeX, sa::GD->sizeY, thickness);
        
        auto layerLV = new G4LogicalVolume(layerS, air, name);

        new G4PVPlacement(nullptr,  // no rotation
            G4ThreeVector(0., 0., displacement),  // translation
            layerLV,  // its logical volume
            name,  // its name
            worldLV,  // its mother volume
            false,  // no boolean operation
            0,  // copy number
            checkOverlaps);  // checking overlaps

        // prepare the replica volumes
        const auto & matName = std::get<0>(param);
        auto * material = LUTmat[matName];

        name = name + std::string("Rep");
        auto layerRepS = new G4Box(name,
            sa::GD->sizeX, sa::GD->sizeY, sa::GD->resZ);
        
        auto layerRepLV = new G4LogicalVolume(layerRepS, material, name);

        int nReplicas = int(thickness / sa::GD->resZ);

        new G4PVReplica(
            name,  // its name
            layerRepLV,  // its logical volume
            layerLV,  // its mother volume
            kZAxis,  // axis of replication
            nReplicas,  // number of replica
            sa::GD->resZ * 2);  // with of replica

        this->logicals.push_back(std::make_pair(name, layerRepLV));
    }
#endif
    return worldPV;
}


void sa::DetectorConstruction::ConstructSDandField()
{
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    for (int i=0; i<sa::GD->layers.size(); i++)
    {   
        std::string SDname = std::string("SD") + std::to_string(i+1);
        auto senseDet = new G4MultiFunctionalDetector(SDname);
        G4SDManager::GetSDMpointer()->AddNewDetector(senseDet);

        G4VPrimitiveScorer* primitive;
        primitive = new G4PSEnergyDeposit("Edep");
        senseDet->RegisterPrimitive(primitive);

        SetSensitiveDetector(std::get<1>(this->logicals[i]), senseDet);
    }
}