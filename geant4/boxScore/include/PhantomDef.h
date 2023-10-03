#ifndef PhantomDef_h
#define PhantomDef_h 1

#include <vector>
#include <string>

namespace bs
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
        std::vector<std::tuple<std::string, double>> layers;
    };

    extern GeomDef* GD;
}

#endif