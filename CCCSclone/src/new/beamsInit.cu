#include <boost/filesystem.hpp>
#include <random>
#include <iostream>

#include "dosecalc_defs.h"
#include "paths.h"
#include "argparse.h"

#include "benchmark_host.cuh"
#include "beam.h"

using namespace old;
namespace fs = boost::filesystem;

float inline dev::rand01() {
    // This function generates a random number between 0 and 1
    return float(std::rand()) / RAND_MAX;
}

int dev::beamsInit(std::vector<BEAM>& beam_arr, int verbose) {
    // Firstly, read a template beam.
    BEAM example;
    std::string beamfile = (fs::path(Paths::Instance()->temp_dir()) / 
        fs::path("omni_beam_list.txt")).string();
    FILE *beam_in;
    if ( (beam_in = fopen(beamfile.c_str(),"r")) == NULL ) {
        printf("Cannot open beam list for reading!\n");
        return -1;
    }

    int beam_check;
    fscanf(beam_in,"%d", &beam_check );

    int beam_id;
    char iso_type[25];
    char iso_loc[25];
    char orient_type[25];
    char fmap_fname[1024];
    float deg_azimuth = 0.0f;
    float deg_zenith = 0.0f;
    float deg_coll = 0.0f;
    fscanf(beam_in,"\n%d %f %f %f %f %s %s %f %f %f %f %f %f %s %f %f %f %f %f %d %d %s",
                    &beam_id, &deg_azimuth, &deg_zenith, &deg_coll, &example.sad,
                    iso_type, iso_loc, &example.isocenter.x, &example.isocenter.y, &example.isocenter.z,
                    &example.source.x, &example.source.y, &example.source.z,
                    orient_type, &example.direction.x, &example.direction.y, &example.direction.z,
                    &example.beamlet_size.x, &example.beamlet_size.y,
                    &example.fmap_size.x, &example.fmap_size.y,
                    fmap_fname );

    example.set_isocenter_type(iso_type);
    example.set_isocenter_location(iso_loc);
    example.set_orientation_type(orient_type);
    example.fmap_fname = std::string(fmap_fname);
    example.azimuth = (deg_azimuth * PI)/180.0f;
    example.zenith = (deg_zenith * PI)/180.0f;
    example.coll = (deg_coll * PI)/180.0f;
    example.uid = 0;

    if (verbose>=2) {
        printf("  Beam Angles   [deg]: (gantry: %3.1f, couch: %3.1f, coll: %3.1f)\n",example.deg_azimuth(),example.deg_zenith(),example.deg_coll());
        printf("  SAD            [cm]: %g\n", example.sad);
        printf("  Beam source    [cm]: (%3.1f, %3.1f, %3.1f)\n",example.source.x,example.source.y,example.source.z);
        printf("  Beam isocenter [cm]: (%3.1f, %3.1f, %3.1f)[%s/%s]\n",example.isocenter.x,example.isocenter.y,example.isocenter.z,example.get_isocenter_type().c_str(), example.get_isocenter_location().c_str());
        printf("  Beam direction:      (%3.1f, %3.1f, %3.1f)[%s]\n",example.direction.x,example.direction.y,example.direction.z,example.get_orientation_type().c_str());
        printf("  Fluence file:        %s\n", example.fmap_fname.c_str() );
        printf("  Beamlet size   [cm]: %g x %g\n",example.beamlet_size.x,example.beamlet_size.y);
        printf("  Fluence dim:         %d x %d\n",example.fmap_size.x,example.fmap_size.y);
    }

    int num_beams = dev::getarg<int>("beamCount");
    beam_arr.clear();
    float zenith_scale = PI;
    float azimuth_scale = PI / 6;
    for (int i=0; i<num_beams; i++) {
        // generate random angles
        beam_arr.emplace_back(example);
        beam_arr.back().zenith = (rand01() - 0.5) * 2 * zenith_scale;
        beam_arr.back().azimuth = (rand01() - 0.5) * 2 * azimuth_scale;
        beam_arr.back().fmap_size = uint2{10, 10};

        if (verbose>=1) {
            std::cout << "beam idx: " << i << ", (zenith, azimuth)=(" << 
                beam_arr.back().deg_zenith() << ", " << beam_arr.back().deg_azimuth() << 
                ") [deg], fluence map: (" << beam_arr.back().fmap_size.x << ", " << 
                beam_arr.back().fmap_size.y << ")" << std::endl;
        }
    }

    return 0;
}