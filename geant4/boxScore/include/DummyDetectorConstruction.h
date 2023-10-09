#ifndef DummyDetectorConstruction_h
#define DummyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4PhysicalVolume;
class G4LogicalVolume;

namespace bs
{
    class DummyDetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        DummyDetectorConstruction() = default;
        ~DummyDetectorConstruction() override = default;

        G4VPhysicalVolume* Construct() override;

        G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    protected:
        G4LogicalVolume* fScoringVolume = nullptr;
    };
}

#endif