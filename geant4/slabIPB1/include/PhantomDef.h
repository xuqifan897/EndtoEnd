#ifndef siPhantomDef_h
#define siPhantomDef_h 1
#include "globals.hh"
#include <vector>

namespace si
{
    class GeomDef
    {
    public:
        GeomDef();
        ~GeomDef();
        GeomDef(GeomDef& old) = delete;
        GeomDef(GeomDef&& old) = delete;
        GeomDef& operator=(const GeomDef right) = delete;

        // display the geometric parameters
        void display();

        // layers are specified by material name and thickness
        std::vector<std::tuple<G4String, G4double, G4double>> layers;
        G4float sizeX;  // half size along X
        G4float sizeY;  // half size along Y
        G4float sizeZ;  // half size along Z, not specified by user, but determined by layers
        G4float resX;  // resolution along X
        G4float resY;  // resolution along Y
        G4float resZ;  // resolution along Z
        G4int dimX;  // dimension along X, not specified by user
        G4int dimY;  // dimension along Y, not specified by user
        G4int dimZ;  // dimension along Z, not specified by user
        G4bool oddXY;  // whether to round dimension to odd integers
    };

    extern GeomDef* GD;
}
#endif