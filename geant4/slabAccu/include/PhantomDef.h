#ifndef PhantomDef_h
#define PhantomDef_h 1

#include <vector>
#include <string>

namespace sa
{
    class GeomDef
    {
    public:
        GeomDef();
        ~GeomDef()=default;
        GeomDef(GeomDef& old) = delete;
        GeomDef(GeomDef&& old) = delete;
        GeomDef& operator=(const GeomDef& right) = delete;

        // display the geometric parameters
        void display();

        // layers are specified by material name and thickness
        std::vector<std::tuple<std::string, double, double>> layers;
        float sizeX;  // half size along X
        float sizeY;  // half size along Y
        float sizeZ;  // half size along Z, not specified by user, but determined by layers
        float resX;  // resolution along X
        float resY;  // resolution along Y
        float resZ;  // resolution along Z
        int dimX;  // dimension along X, not specified by user
        int dimY;  // dimension along Y, not specified by user
        int dimZ;  // dimension along Z, not specified by user
        bool oddXY;  // whether to round dimension to odd integers
    };

    extern GeomDef* GD;
}

#endif