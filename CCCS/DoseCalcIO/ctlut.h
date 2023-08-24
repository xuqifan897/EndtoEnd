#ifndef __CTLUT_H__
#define __CTLUT_H__

#include <iostream>
#include <string>
#include <vector>
#include <utility>


// represents a single point in the lookup table
struct LUTPOINT {
    LUTPOINT(std::string label, float hunits, float massdens, float reledens=0)
        : label(label), hunits(hunits), massdens(massdens), reledens(reledens) {}

    std::string label;
    float hunits;     // [HU]
    float massdens;   // [g/cm^3]
    float reledens;   // electron density relative to water
};

// collection of LUTPOINTS with convenience methods and interpolation
struct CTLUT {
    enum class INTERPTYPE {
        LINEAR,
    };

    CTLUT(INTERPTYPE interp_style=INTERPTYPE::LINEAR) : label(""), interp_style(interp_style) {};
    CTLUT(std::string label, INTERPTYPE interp_style=INTERPTYPE::LINEAR)
        : label(label), interp_style(interp_style) {}

    std::string label;
    INTERPTYPE interp_style;

    std::vector<LUTPOINT> points;
    void sort();

    friend std::ostream& operator<<(std::ostream& os, const CTLUT& ctlut);
};

int load_lookup_table(CTLUT& lut, std::string filepath, int verbose=0);

#endif // __CTLUT_H__

