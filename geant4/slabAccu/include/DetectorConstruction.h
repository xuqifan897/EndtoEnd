#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace sa
{
    class DetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DetectorConstruction() = default;
        ~DetectorConstruction() override = default;

        G4VPhysicalVolume* Construct() override;
        void ConstructSDandField() override;
    
    private:
        std::vector<std::pair<std::string, G4LogicalVolume*>> logicals;
        G4LogicalVolume* logicWorld;
    };
}

#endif