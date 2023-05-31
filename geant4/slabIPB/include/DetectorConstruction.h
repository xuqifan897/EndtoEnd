#ifndef siDetectorConstruction_h
#define siDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

namespace si
{
    class DetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DetectorConstruction();
        ~DetectorConstruction();
        DetectorConstruction(DetectorConstruction& old) = delete;
        DetectorConstruction(DetectorConstruction&& old) = delete;
        DetectorConstruction& operator=(const DetectorConstruction old) = delete;

        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();
    
    private:
        // name of logical volume, logical volume pointer
        std::vector<std::pair<G4String, G4LogicalVolume*>> logicals;
    };
}

#endif