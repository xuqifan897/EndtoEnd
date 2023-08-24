#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "dosecalc-preprocess"

#include <boost/test/unit_test.hpp>
#include "../dosecalc-preprocess/fmapProcessing.hpp"
#include "helper_math.h"

#include <cstdio>
BOOST_AUTO_TEST_CASE(constructors) {
    std::vector<float> fmap = {1,0,1,
                               0,1,0,
                               1,1,0};
    uint2 fmap_size = uint2{3,3};
    BOOST_CHECK(!fmap_is_apertureready(fmap,fmap_size));
    fmap_post_apertureready(fmap, fmap_size);
    BOOST_CHECK(fmap_is_apertureready(fmap,fmap_size));
    std::vector<float> expected = {1,1,1,
                                   0,1,0,
                                   1,1,0};
    BOOST_CHECK_EQUAL_COLLECTIONS(fmap.begin(), fmap.end(), expected.begin(), expected.end());

    fmap = {1,0,1,0,1,
            1,0,0,1,0,
            0,0,1,0,1};
    fmap_size = uint2{5,5};
    BOOST_CHECK(!fmap_is_apertureready(fmap,fmap_size));
    fmap_post_apertureready(fmap, fmap_size);
    BOOST_CHECK(fmap_is_apertureready(fmap,fmap_size));
    expected = {1,1,1,1,1,
                1,1,1,1,0,
                0,0,1,1,1};
    BOOST_CHECK_EQUAL_COLLECTIONS(fmap.begin(), fmap.end(), expected.begin(), expected.end());
}
