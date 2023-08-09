#include <vector>
#include <string>

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "PhantomDef.h"
#include "config.h"

ts::GeomDef* ts::GD = nullptr;

ts::GeomDef::GeomDef()
{
    this->layers_nominal = std::vector<std::tuple<std::string, double>>();

    // please specify the parameters in the below block
//000000000000000000000000000000000000000000000000000000000000000000000000
# if PHANTOM == 0
    this->layers_nominal.push_back(std::make_tuple("adipose", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("lung", 4.8*cm));
    this->layers_nominal.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("adipose", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("bone", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("muscle", 0.8*cm));
    this->layers_nominal.push_back(std::make_tuple("adipose", 0.8*cm));
#endif

    this->radius = 10. * cm;
    this->resR = 0.05 * cm;
    this->resZ = 0.05 * cm;
//00000000000000000000000000000000000000000000000000000000000000000000000

    this->dimR = std::round(this->radius / this->resR);
    this->radius = this->dimR * this->resR;

    // firstly, we ensure that the thickness of each layer is a multiple of resZ
    this->sizeZ = 0;
    for (auto it=this->layers_nominal.begin(); 
        it!=this->layers_nominal.end(); it++)
    {
        float thickness = std::get<1>(*it);
        thickness = std::round(thickness / this->resZ) * this->resZ;
        std::get<1>(*it) = thickness;
        this->sizeZ += thickness;
    }

    // Then, we construct physical layers using the information of the nominal layers
    this->dimZ = std::round(this->sizeZ / this->resZ);
    this->layers_physical = std::vector<std::tuple<
        std::string, double, double>>(this->dimZ);
    
    int i=0;
    for (auto it=this->layers_physical.begin(); 
        it!=this->layers_physical.end(); it++)
    {
        double slabThickness = std::get<1>(*it);
        std::string& slabMaterial = std::get<0>(*it);
        int slabDim = std::round(slabThickness / this->resZ);
        for (int j=0; j<slabDim; j++)
        {
            std::get<0>(this->layers_physical[i+j]) = slabMaterial;
            std::get<1>(this->layers_physical[i+j]) = this->resZ;
            std::get<2>(this->layers_physical[i+j]) = 
                (2*(i+j)+1) * this->resZ - this->sizeZ;
        }
        i += slabDim;
    }

    display();
}


void ts::GeomDef::display()
{   
    int count = 1;
    std::cout << "Nominal layers:" << std::endl;
    for (auto it=this->layers_nominal.begin(); 
        it!=this->layers_nominal.end(); it++)
    {
        std::cout << "layer: " << count << ", material: " << std::get<0>(*it);
    }
}