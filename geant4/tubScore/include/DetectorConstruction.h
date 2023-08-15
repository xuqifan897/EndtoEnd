#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace ts
{
    class DetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DetectorConstruction() = default;
        ~DetectorConstruction() = default;

        G4VPhysicalVolume* Construct() override;
        void ConstructSDandField() override;
    
    private:
        std::vector<G4LogicalVolume*> logicals;
    };
}

#endif