#include "G4SystemOfUnits.hh"
#include <vector>
#include "PhantomDef.h"
#include "config.h"

si::GeomDef* si::GD = nullptr;

si::GeomDef::GeomDef()
{
    // please specify the parameters below
// 000000000000000000000000000000000000000000000000000000000000000000000000000
    // It is noted that, all numbers are half thicknesses
    // In the format of (material, thickness, offset)
    this->layers = std::vector<std::tuple<G4String, G4double, G4double>>();
#if PHANTOM == SLAB
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("lung", 4.8*cm, 0.));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm, 0.));
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm, 0.));
#endif

# if PHANTOM == WATER
    this->layers.push_back(std::make_tuple("water", 12.8*cm, 0.));
#endif

    // Half size
    this->sizeX = 10.0 * cm;
    this->sizeY = 10.0 * cm;

    // Half resolution
    this->resX = 0.05 * cm;
    this->resY = 0.05 * cm;
    this->resZ = 0.05 * cm;

    // whether to use odd dimensions along X and Y dimension
    this->oddXY = true;
// 00000000000000000000000000000000000000000000000000000000000000000000000000000

    // The folloowing parameters should be calculated instead of being specified
    this->sizeZ = 0.;
    for (auto it = this->layers.begin(); it != this->layers.end(); it++)
        this->sizeZ += std::get<1>(*it);
    
    // calculate offset
    G4double thickCum = 0;
    for (auto it = this->layers.begin(); it != this->layers.end(); it++)
    {
        std::get<2>(*it) = thickCum + std::get<1>(*it) - this->sizeZ;
        thickCum += 2 * std::get<1>(*it);
    }

    this->dimX = int(this->sizeX / this->resX);
    this->dimY = int(this->sizeY / this->resY);
    this->dimZ = int(this->sizeZ / this->resZ);
    if (this->oddXY)
    {
        if (this->dimX % 2 == 0)
            this->dimX -= 1;
        if (this->dimY % 2 == 0)
            this->dimY -= 1;
    }
    else
    {
        if (this->dimX % 2 == 1)
            this->dimX -= 1;
        if (this->dimY % 2 == 1)
            this->dimY -= 1;
    }

    this->sizeX = this->dimX * this->resX;
    this->sizeY = this->dimY * this->resY;
    this->sizeZ = this->dimZ * this->resZ;

    display();
}

void si::GeomDef::display()
{
    G4cout << "Phantom geometry:" << G4endl;
    G4cout << "Size: (" << std::setprecision(4) << this->sizeX / cm << ", " 
        << std::setprecision(4) << this->sizeY / cm << ", " 
        << std::setprecision(4) << this->sizeZ / cm << ")[cm]" << G4endl;
    G4cout << "Dimension: (" << this->dimX << ", " << this->dimY << ", "
        << this->dimZ << ")" << G4endl;
    G4cout << "Resolution: (" << this->resX / cm << ", " << this->resY / cm 
        << ", " << this->resZ / cm << ")[cm]" << G4endl;
    G4cout << "Layers:" << G4endl;
    G4cout << std::setw(10) << std::right << "material" 
        << std::setw(10) << std::right << "thickness" 
        << std::setw(10) << std::right << "offset" << G4endl;
    for (auto it = this->layers.begin(); it != this->layers.end(); it++)
    {
        G4cout << std::setw(10) << std::right << std::get<0>(*it) 
            << std::setw(10) << std::right << std::setprecision(3) 
            << std::get<1>(*it)/cm << "[cm]"
            << std::setw(10) << std::right << std::setprecision(3)
            << std::get<2>(*it)/cm << "[cm]" << G4endl;
    }
}