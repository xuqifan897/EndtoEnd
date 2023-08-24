#ifndef PhantomDef_h
#define PhantomDef_h 1

#include <vector>
#include <string>

namespace ts
{
    class GeomDef
    {
    public:
        GeomDef();
        ~GeomDef() = default;
        GeomDef(GeomDef& old) = delete;
        GeomDef(GeomDef&& old) = delete;
        GeomDef& operator=(const GeomDef& right) = delete;

        // display the geometric parameters
        void display();

        // layers are specified by material name and thickness
        // (material, half thickness)
        std::vector<std::tuple<std::string, double>> layers_nominal;
        // (material, half thickness, offset)
        std::vector<std::tuple<std::string, double, double>> layers_physical;

        float radius;  // The radius of the disc slab
        float sizeZ; // The half dimension along z

        float resR; // The resolution of the radius
        float resZ; // The resolution along z

        int dimR; // The number of replicas in the radius dimension
        int dimZ; // The number of replicas in the z direction
    };

    extern GeomDef* GD;
}

#endif