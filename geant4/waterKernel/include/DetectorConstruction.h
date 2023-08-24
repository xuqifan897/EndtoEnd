#ifndef wkDetectorConstruction_h
#define wkDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

namespace wk
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
        float sizeX;
        float sizeY;
        float sizeZ;

        G4LogicalVolume* experimentalHall_log;
    };
}

#endif