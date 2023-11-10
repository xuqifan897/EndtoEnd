#include "beam.h"
#include "geometry.h"
#include "math.h"
#include "paths.h"

#include <string>
#include <memory>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cmath>
#include <helper_cuda.h>
#include <helper_math.h>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "dosecalc_defs.h"

float DEFAULT_SAD = 100.f;

void old::BEAM::set_isocenter_type(const std::string& iso_type) {
    if (iso_type == "man")
        isocenter_type = ISO_T::MANUAL;
    else if (iso_type == "ptv")
        isocenter_type = ISO_T::PTV_CENTROID;
    else
        isocenter_type = ISO_T::UNSPEC;
}
std::string old::BEAM::get_isocenter_type() const {
    std::string iso_type;
    switch (isocenter_type) {
        case BEAM::ISO_T::UNSPEC :
            iso_type="unspec"; break;
        case BEAM::ISO_T::PTV_CENTROID :
            iso_type="ptv"; break;
        case BEAM::ISO_T::MANUAL :
            iso_type="man"; break;
        default : iso_type="unknown"; break;
    }
    return iso_type;
}
void old::BEAM::set_isocenter_location(const std::string& iso_loc) {
    if (iso_loc == "in")
        isocenter_location = ISO_LOC_T::IN_PTV;
    else if (iso_loc == "out")
        isocenter_location = ISO_LOC_T::OUT_PTV;
    else
        isocenter_location = ISO_LOC_T::UNSPEC;
}
std::string old::BEAM::get_isocenter_location() const {
    std::string iso_loc;
    switch (isocenter_location) {
        case BEAM::ISO_LOC_T::UNSPEC :
            iso_loc="unspec"; break;
        case BEAM::ISO_LOC_T::IN_PTV :
            iso_loc="in"; break;
        case BEAM::ISO_LOC_T::OUT_PTV :
            iso_loc="out"; break;
        default : iso_loc="unknown"; break;
    }
    return iso_loc;
}
void old::BEAM::set_orientation_type(const std::string& orient_type) {
    if (orient_type == "man" || orient_type=="dir" || orient_type=="direction")
        orientation_type = ORIENT_T::DIRECTION;
    else if (orient_type == "auto" || orient_type=="angle")
        orientation_type = ORIENT_T::ANGLE;
    else
        orientation_type = ORIENT_T::UNSPEC;
}
std::string old::BEAM::get_orientation_type() const {
    std::string orient_type;
    switch (orientation_type) {
        case ORIENT_T::DIRECTION :
            orient_type="man"; break;
        case ORIENT_T::ANGLE :
            orient_type="auto"; break;
        default : orient_type="unknown"; break;
    }
    return orient_type;
}
#define eps 1e-6
CUDEV_FXN float3 old::BEAM::calc_source_from_angles(float gantry_rot_rad, float couch_rot_rad, float3 iso, float sad) {
    // float3 src = iso + make_float3( 0.f, -sad, 0.f );
    // return inverseRotateBeamRHS( src, iso, gantry_rot_rad, couch_rot_rad, 0.f);
    float3 src = make_float3(0.f, -sad, 0.f);
    return inverseRotateBeamAtOriginRHS(src, gantry_rot_rad, couch_rot_rad, 0.f) + iso;
}
CUDEV_FXN float2 old::BEAM::calc_angles_from_dir(float3 dir, float3 iso, float sad) {
    // Calculate gantry+couch angles for a given beam direction vector
    float3 ndir = -1*normalize(dir);
    float theta = acosf(-ndir.y);
    float phi = atanf(ndir.z/ndir.x);

    // resolve ambiguous directions (interface of quadrants)
    if (closeto(ndir.z, 0.f)) { phi = 0.f; }
    else if (closeto(ndir.x, 0.f)) { phi = PI/2.f; }

    // make gantry rotation consistent with couch angles
    if (closeto(ndir.x, 0.f) && closeto(ndir.z, 0.f)) { theta = (ndir.y<=(-eps)?0:PI); }
    else if (ndir.x<(-eps) || (closeto(ndir.x, 0.f) && ndir.z < (-eps))) { theta *= -1; }

    return float2{ theta, phi };
}

int old::load_omni_beam_list(std::vector<BEAM>& beams, int beam_count, int verbose)
{
    std::string beamfile = (fs::path(Paths::Instance()->temp_dir()) / 
        fs::path("omni_beam_list.txt")).string();
    FILE *beam_in;
    if ( (beam_in = fopen(beamfile.c_str(),"r")) == NULL ) {
        printf("Cannot open beam list for reading!\n");
        return -1;
    }
    int beam_check;
    fscanf(beam_in,"%d", &beam_check );
    if (beam_count < beam_check) { beam_check = beam_count; };
    if (verbose>=2) {
        printf("BEAM-DATA:\n");
    }
    for (int b=0; b<beam_check; b++) {
        int beam_id;
        char iso_type[25];
        char iso_loc[25];
        char orient_type[25];
        char fmap_fname[1024];
        float deg_azimuth = 0.0f;
        float deg_zenith = 0.0f;
        float deg_coll = 0.0f;
        fscanf(beam_in,"\n%d %f %f %f %f %s %s %f %f %f %f %f %f %s %f %f %f %f %f %d %d %s",
                        &beam_id, &deg_azimuth, &deg_zenith, &deg_coll, &beams[b].sad,
                        iso_type, iso_loc, &beams[b].isocenter.x, &beams[b].isocenter.y, &beams[b].isocenter.z,
                        &beams[b].source.x, &beams[b].source.y, &beams[b].source.z,
                        orient_type, &beams[b].direction.x, &beams[b].direction.y, &beams[b].direction.z,
                        &beams[b].beamlet_size.x, &beams[b].beamlet_size.y,
                        &beams[b].fmap_size.x, &beams[b].fmap_size.y,
                        fmap_fname );

        beams[b].set_isocenter_type(iso_type);
        beams[b].set_isocenter_location(iso_loc);
        beams[b].set_orientation_type(orient_type);
        beams[b].fmap_fname = std::string(fmap_fname);
        beams[b].azimuth = (deg_azimuth * PI)/180.0f;
        beams[b].zenith = (deg_zenith * PI)/180.0f;
        beams[b].coll = (deg_coll * PI)/180.0f;
        beams[b].uid = b;

        if (verbose>=2) {
            printf("  Beam Angles   [deg]: (gantry: %3.1f, couch: %3.1f, coll: %3.1f)\n",beams[b].deg_azimuth(),beams[b].deg_zenith(),beams[b].deg_coll());
            printf("  SAD            [cm]: %g\n", beams[b].sad);
            printf("  Beam source    [cm]: (%3.1f, %3.1f, %3.1f)\n",beams[b].source.x,beams[b].source.y,beams[b].source.z);
            printf("  Beam isocenter [cm]: (%3.1f, %3.1f, %3.1f)[%s/%s]\n",beams[b].isocenter.x,beams[b].isocenter.y,beams[b].isocenter.z,beams[b].get_isocenter_type().c_str(), beams[b].get_isocenter_location().c_str());
            printf("  Beam direction:      (%3.1f, %3.1f, %3.1f)[%s]\n",beams[b].direction.x,beams[b].direction.y,beams[b].direction.z,beams[b].get_orientation_type().c_str());
            printf("  Fluence file:        %s\n", beams[b].fmap_fname.c_str() );
            printf("  Beamlet size   [cm]: %g x %g\n",beams[b].beamlet_size.x,beams[b].beamlet_size.y);
            printf("  Fluence dim:         %d x %d\n",beams[b].fmap_size.x,beams[b].fmap_size.y);
        }

        // // enable all beamlets, no need to allocate fluence map
        // beams[b].fluence_map = std::vector<float>(beams[b].fmap_size.x * beams[b].fmap_size.y);
        // if (verbose) { printf("  "); }
        // load_fluence_map( beams[b].fluence_map, beams[b].fmap_size, beams[b].fmap_fname );
        // if (verbose) { printf("\n"); }
    }
    fclose(beam_in);
    return 0;
}