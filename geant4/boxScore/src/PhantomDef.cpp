#include <vector>
#include <string>

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "PhantomDef.h"

bs::GeomDef* bs::GD = nullptr;

bs::GeomDef::GeomDef()
{
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers.push_back(std::make_tuple("lung", 4.8*cm));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm));
    this->layers.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers.push_back(std::make_tuple("adipose", 0.8*cm));
}

void bs::GeomDef::display()
{
    for (int i=0; i<this->layers.size(); i++)
    {
        std::string& mat = std::get<0>(this->layers[i]);
        float thickness = std::get<1>(this->layers[i]);
        std::cout << "material:" << std::setw(15) << std::right << mat << ", ";
        std::cout << "thickness:" << std::setw(14) << std::right << thickness/cm << "cm" << std::endl;
    }
}