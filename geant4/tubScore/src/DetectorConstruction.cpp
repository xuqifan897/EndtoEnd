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

G4VPhysicalVolume* ts::DetectorConstruction::Construct()
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

    // display material information
    std::cout << std::endl;
    std::cout << "Material information: " << std::endl;

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
    auto worldS = new G4Tubs("World", 0., GD->radius, 
        GD->sizeZ, 0.*deg, 360.*deg);  // size
    
    auto worldLV = new G4LogicalVolume(
        worldS,  // solid
        air,  // material
        "World");  // name
    
    auto worldPV = new G4PVPlacement(
        nullptr,  // no rotation
        G4ThreeVector(),  // at (0, 0, 0)
        worldLV,  // its logical volume
        "World",  // its name
        nullptr,  // no mother volume
        false,  // no boolean operator
        0,  // copy number
        checkOverlaps);  // checking overlaps

    int count = 1;
    for (auto it=GD->layers_physical.begin(); 
        it!=GD->layers_physical.end(); it++)
    {
        const std::string & material = std::get<0>(*it);
        const double & thickness = std::get<1>(*it);
        const double & displacement = std::get<2>(*it);
        auto* mat = LUTmat[material];

        std::string name = "layer" + std::to_string(count);
        count ++;
        auto layerS = new G4Tubs(name, 0., GD->radius, 
            thickness, 0.*deg, 360.*deg);
        auto layerLV = new G4LogicalVolume(layerS, mat, name);
        // auto layerPV =
            new G4PVPlacement(
                nullptr,  // no rotation
                G4ThreeVector(0., 0., displacement),  // displacement
                layerLV,  // its logical volume
                name,  // its name
                worldLV,  // its mother volume
                false,  // no boolean operation
                0,  // copy number
                checkOverlaps);  // checking overlaps
        
        // prepare the replica volumes
        name = name + std::string("Rep");
        auto layerRepS = new G4Tubs(name, 0., GD->radius, 
            thickness, 0.*deg, 360.*deg);
        auto layerResLV = new G4LogicalVolume(layerRepS, mat, name);
        // auto layerResPV = 
            new G4PVReplica(name, layerResLV, layerLV, 
                kRho, GD->dimR, GD->resR);
        
        this->logicals.push_back(layerResLV);
    }
    return worldPV;
}

void ts::DetectorConstruction::ConstructSDandField()
{
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    for (int i=0; i<this->logicals.size(); i++)
    {
        std::string SDname = std::string("SD") + std::to_string(i+1);
        auto senseDet = new G4MultiFunctionalDetector(SDname);
        G4SDManager::GetSDMpointer()->AddNewDetector(senseDet);

        G4VPrimitiveScorer* primitive;
        primitive = new G4PSEnergyDeposit("Edep");
        senseDet->RegisterPrimitive(primitive);

        SetSensitiveDetector(this->logicals[i], senseDet);
    }
}