#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace wp
{
    class DetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DetectorConstruction() = default;
        ~DetectorConstruction() = default;

        G4VPhysicalVolume* Construct() override;
        void ConstructSDandField() override;
    
    private:
        // The smallest detector element
        // G4LogicalVolume* detElement;
        std::vector<G4LogicalVolume*> SenseDetList;
    };
}

#endif