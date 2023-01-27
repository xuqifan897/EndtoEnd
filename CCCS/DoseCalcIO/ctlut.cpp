#include "ctlut.h"

#include <cstdio>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <algorithm>

#include "./io_helpers.h"

using namespace dcio;

std::ostream& operator<<(std::ostream& os, const CTLUT& ctlut) {
    os << "CT Lookup Table ("<<ctlut.label<<"):"<<std::endl;
    os << "---------------------------------------------" << std::endl;
    if (!ctlut.points.size()) {
        os << "  EMPTY" << std::endl;
    } else {
        os << "   #     HU#    g/cm3   rel   label          " << std::endl;
        os << "  ---  -------  -----  -----  ---------------" << std::endl;
        char buf[100];
        int ii = 0;
        for (const auto& pt : ctlut.points) {
            ii++;
            sprintf(buf, "  %3d  %7.1f  %5.3f  %5.3f  %s", ii, pt.hunits, pt.massdens, pt.reledens, pt.label.c_str());
            os << buf << std::endl;
        }
    }
    return os;
}

void CTLUT::sort() {
    std::sort(points.begin(), points.end(),
            [] (const LUTPOINT& p1, const LUTPOINT& p2) { return (p1.hunits < p2.hunits); }
            );
}


int load_lookup_table(CTLUT& lut, std::string filepath, int verbose) {
    /* LUT file format:
     * Each line should contain the material name, CT number, corresponding electron density, and optionally the density relative to water (currently unused)
     *    name   CT#      density  rel_density
     *    <str>  <float>  <float>  <float>
     * lines beginning with "#" are ignored as comments, empty lines between entries are ignored, lines are read until EOF
     * */
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Cannot open lookup table file: \""<<filepath<<"\" for reading" << std::endl;
        std::cout << "failed with error ("<<errno<<"): " << std::strerror(errno) << std::endl;
        return false;
    }

    std::string line;
    int currentline = 0;
    while (std::getline(file, line)) {
        currentline++;

        // test for empty
        if (line.empty()) {
            if (verbose>=2) { std::cout << "ignoring empty line: "<<currentline<<std::endl; }
            continue;
        }

        // test for comment
        if (is_comment_string(line, '#')) {
            if (verbose>=2) { std::cout << "ignoring comment on line "<<currentline<<std::endl; }
            continue;
        }

        // split line into fields
        std::vector<std::string> fields{};
        tokenize_string(line, fields, " ");
        if (fields.size() < 3 ) {
            std::cout << "LUT spec on line " << currentline << " is invalid. Please check documentation for valid specification" << std::endl;
            return false;
        }

        // parse fields
        lut.points.emplace_back(fields[0], std::stof(fields[1]), std::stof(fields[2]));
        if (verbose>=2) {
            std::cout << "added LUT point on line " << currentline<< ": "<<
                lut.points.back().label<<" "<<lut.points.back().hunits<<" "<<lut.points.back().massdens<< std::endl;
        }
    }

    return true;
}
