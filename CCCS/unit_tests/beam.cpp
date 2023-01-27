#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "dosecalcIO_beam"

#include <boost/test/unit_test.hpp>
#include "../DoseCalcIO/beam.h"
#include "helper_math.h"
#include "../include/dosecalc_defs.h"
#include "../Utilities/math.h"

#include <iostream>
#include <cstdio>
#include <vector>
namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_CASE(constructors) {
    std::vector<float3> angles;
    // coplanar angles
    for (int ii=0; ii<12; ii++) {
        angles.push_back({45.f*(float)ii,0});
    }

    // non-coplanar angles
    for (int jj=0; jj<12; jj++) {
        for (int ii=0; ii<12; ii++) {
            angles.push_back({45.f*(float)ii,45.f*(float)jj});
        }
    }

    for (const float3& ang : angles) {
        printf("testing gantry:%0.1f couch:%0.1f\n", ang.x, ang.y);
        BEAM beam1{};
        beam1.azimuth = ang.x*PI/180.0;
        beam1.zenith = ang.y*PI/180.0;
        beam1.isocenter = float3{0,0,0};
        beam1.sad = 100;
        beam1.reconfigure();
        std::cout << beam1 << std::endl;

        BEAM beam2{};
        beam2.direction = beam1.direction;
        beam2.isocenter = float3{0,0,0};
        beam2.orientation_type = BEAM::ORIENT_T::DIRECTION;
        beam2.sad = 100;
        beam2.reconfigure();
        std::cout << beam2 << std::endl;

        BEAM beam3{};
        beam3.azimuth = beam2.azimuth;
        beam3.zenith  = beam2.zenith;
        beam3.sad = 100;
        beam3.reconfigure();
        std::cout << beam3 << std::endl;

        printf("  dir: %0.2f %0.2f %0.2f\n", beam1.direction.x, beam1.direction.y, beam1.direction.z);
        // BOOST_CHECK_MESSAGE(closeto(beam1.azimuth, beam2.azimuth),
        //         "Gantry angle mismatch (ang:"<<beam1.azimuth<<", dir:"<<beam2.azimuth<<")");
        // BOOST_CHECK_MESSAGE(closeto(beam1.zenith, beam2.zenith),
        //         "Couch angle mismatch (ang:"<<beam1.zenith<<", dir:"<<beam2.zenith<<")");
        BOOST_CHECK_MESSAGE(closeto(beam1.direction.x, beam3.direction.x),
                "Direction (X) mismatch (ang:"<<beam1.direction.x<<", dir:"<<beam3.direction.x<<")");
        BOOST_CHECK_MESSAGE(closeto(beam1.direction.y, beam3.direction.y),
                "Direction (Y) mismatch (ang:"<<beam1.direction.y<<", dir:"<<beam3.direction.y<<")");
        BOOST_CHECK_MESSAGE(closeto(beam1.direction.z, beam3.direction.z),
                "Direction (Z) mismatch (ang:"<<beam1.direction.z<<", dir:"<<beam3.direction.z<<")");
    }
}
