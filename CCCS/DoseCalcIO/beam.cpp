#include "beam.h"

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

#include "dosecalc_defs.h"
#include "Utilities/logging.h"
#include "./io_helpers.h"
#include "./fluence.h"
#include "./paths.h"
#include "Utilities/math.h"
#include "CudaUtilities/geometry.cuh"

using namespace dcio;

float DEFAULT_SAD = 100.f;  // used when sad: %f isn't specified for beam in beamfile [unit: cm]

///////// BEAM METHOD DEFINITIONS //////////
void BEAM::set_isocenter_type(const std::string& iso_type) {
    if (iso_type == "man")
        isocenter_type = ISO_T::MANUAL;
    else if (iso_type == "ptv")
        isocenter_type = ISO_T::PTV_CENTROID;
    else
        isocenter_type = ISO_T::UNSPEC;
}
std::string BEAM::get_isocenter_type() const {
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
void BEAM::set_isocenter_location(const std::string& iso_loc) {
    if (iso_loc == "in")
        isocenter_location = ISO_LOC_T::IN_PTV;
    else if (iso_loc == "out")
        isocenter_location = ISO_LOC_T::OUT_PTV;
    else
        isocenter_location = ISO_LOC_T::UNSPEC;
}
std::string BEAM::get_isocenter_location() const {
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
void BEAM::set_orientation_type(const std::string& orient_type) {
    if (orient_type == "man" || orient_type=="dir" || orient_type=="direction")
        orientation_type = ORIENT_T::DIRECTION;
    else if (orient_type == "auto" || orient_type=="angle")
        orientation_type = ORIENT_T::ANGLE;
    else
        orientation_type = ORIENT_T::UNSPEC;
}
std::string BEAM::get_orientation_type() const {
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
CUDEV_FXN float3 BEAM::calc_source_from_angles(float gantry_rot_rad, float couch_rot_rad, float3 iso, float sad) {
    // float3 src = iso + make_float3( 0.f, -sad, 0.f );
    // return inverseRotateBeamRHS( src, iso, gantry_rot_rad, couch_rot_rad, 0.f);
    float3 src = make_float3(0.f, -sad, 0.f);
    return inverseRotateBeamAtOriginRHS(src, gantry_rot_rad, couch_rot_rad, 0.f) + iso;
}
CUDEV_FXN float2 BEAM::calc_angles_from_dir(float3 dir, float3 iso, float sad) {
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
CUDEV_FXN float3 BEAM::calc_source_from_dir(float3 dir, float3 iso, float sad) {
    return iso - normalize(dir) * sad;
}
CUDEV_FXN float3 BEAM::calc_dir_from_source(float3 iso, float3 source) {
    return normalize(iso-source);
}

CUDEV_FXN float3 BEAM::calc_source_from_angles() const {
    return BEAM::calc_source_from_angles(azimuth, zenith, isocenter, sad);
}
CUDEV_FXN float2 BEAM::calc_angles_from_dir()    const {
    return BEAM::calc_angles_from_dir(direction, isocenter, sad);
}
CUDEV_FXN float3 BEAM::calc_source_from_dir()    const {
    return BEAM::calc_source_from_dir(direction, isocenter, sad);
}
CUDEV_FXN float3 BEAM::calc_dir_from_source()    const {
    return BEAM::calc_dir_from_source(isocenter, source);
}

H5::CompType BEAM::getFileCompoundType() {
    // tuple dataspaces
    hsize_t tuple3_dims[] = { 3 };
    H5::DataSpace tuple3(1, tuple3_dims);
    hsize_t tuple2_dims[] = { 2 };
    H5::DataSpace tuple2(1, tuple2_dims);

    // define file datatype - user preference
    size_t dataset_size = sizeof(COMP_BEAM_T);// + sizeof(float)*(fmap_size.x*fmap_size.y);
    H5::CompType beam_file_t( dataset_size );

    // define elementary dattypes
    H5::ArrayType src_t(H5::PredType::IEEE_F32LE, 1, tuple3_dims);
    H5::ArrayType dir_t(H5::PredType::IEEE_F32LE, 1, tuple3_dims);
    H5::ArrayType iso_t(H5::PredType::IEEE_F32LE, 1, tuple3_dims);
    H5::ArrayType fmap_dims_t(H5::PredType::STD_U32LE, 1, tuple2_dims);
    H5::ArrayType beamlet_size_t(H5::PredType::IEEE_F32LE, 1, tuple2_dims);

    beam_file_t.insertMember("uid",             HOFFSET(COMP_BEAM_T, uid),            H5::PredType::STD_U16LE);
    beam_file_t.insertMember("gantry_rot_rad",  HOFFSET(COMP_BEAM_T, gantry_rot_rad), H5::PredType::IEEE_F32LE);
    beam_file_t.insertMember("couch_rot_rad",   HOFFSET(COMP_BEAM_T, couch_rot_rad),  H5::PredType::IEEE_F32LE);
    beam_file_t.insertMember("coll_rot_rad",    HOFFSET(COMP_BEAM_T, coll_rot_rad),   H5::PredType::IEEE_F32LE);
    beam_file_t.insertMember("src_coords_cm",   HOFFSET(COMP_BEAM_T, src_coords_cm),  src_t);
    beam_file_t.insertMember("direction",       HOFFSET(COMP_BEAM_T, direction),      dir_t);
    beam_file_t.insertMember("iso_coords_cm",   HOFFSET(COMP_BEAM_T, iso_coords_cm),  iso_t);
    beam_file_t.insertMember("fmap_dims",       HOFFSET(COMP_BEAM_T, fmap_dims),      fmap_dims_t);
    beam_file_t.insertMember("beamlet_size_cm", HOFFSET(COMP_BEAM_T, beamlet_size),   beamlet_size_t);

    return beam_file_t;
}
H5::CompType BEAM::getMemCompoundType() {
    // tuple dataspaces
    hsize_t tuple3_dims[] = { 3 };
    H5::DataSpace tuple3(1, tuple3_dims);
    hsize_t tuple2_dims[] = { 2 };
    H5::DataSpace tuple2(1, tuple2_dims);

    // Define memory datatype - necessary to read from mem buffer properly when writing
    size_t dataset_size = sizeof(COMP_BEAM_T);// + sizeof(float)*(fmap_size.x*fmap_size.y);
    H5::CompType beam_mem_t( dataset_size );

    // define elementary dattypes
    H5::ArrayType src_t(H5::PredType::NATIVE_FLOAT, 1, tuple3_dims);
    H5::ArrayType dir_t(H5::PredType::NATIVE_FLOAT, 1, tuple3_dims);
    H5::ArrayType iso_t(H5::PredType::NATIVE_FLOAT, 1, tuple3_dims);
    H5::ArrayType fmap_dims_t(H5::PredType::NATIVE_UINT, 1, tuple2_dims);
    H5::ArrayType beamlet_size_t(H5::PredType::NATIVE_FLOAT, 1, tuple2_dims);

    beam_mem_t.insertMember("uid",             HOFFSET(COMP_BEAM_T, uid),            H5::PredType::NATIVE_USHORT);
    beam_mem_t.insertMember("gantry_rot_rad",  HOFFSET(COMP_BEAM_T, gantry_rot_rad), H5::PredType::NATIVE_FLOAT);
    beam_mem_t.insertMember("couch_rot_rad",   HOFFSET(COMP_BEAM_T, couch_rot_rad),  H5::PredType::NATIVE_FLOAT);
    beam_mem_t.insertMember("coll_rot_rad",    HOFFSET(COMP_BEAM_T, coll_rot_rad),   H5::PredType::NATIVE_FLOAT);
    beam_mem_t.insertMember("src_coords_cm",   HOFFSET(COMP_BEAM_T, src_coords_cm),  src_t);
    beam_mem_t.insertMember("direction",       HOFFSET(COMP_BEAM_T, direction),      dir_t);
    beam_mem_t.insertMember("iso_coords_cm",   HOFFSET(COMP_BEAM_T, iso_coords_cm),  iso_t);
    beam_mem_t.insertMember("fmap_dims",       HOFFSET(COMP_BEAM_T, fmap_dims),      fmap_dims_t);
    beam_mem_t.insertMember("beamlet_size_cm", HOFFSET(COMP_BEAM_T, beamlet_size),   beamlet_size_t);

    return beam_mem_t;
}
int BEAM::_readFromHDF5(BEAM& beam, H5::Group& h5group) {

    // get compound representations
    H5::CompType beam_mem_t = BEAM::getMemCompoundType();

    // read compound attributes
    {
        COMP_BEAM_T beam_transfer {};
        auto attr = h5group.openAttribute("beam_specs");
        attr.read(beam_mem_t, &beam_transfer);

        beam.uid     = beam_transfer.uid;
        beam.azimuth = beam_transfer.gantry_rot_rad;
        beam.zenith  = beam_transfer.couch_rot_rad;
        beam.coll    = beam_transfer.coll_rot_rad;
        ARR3VECT(beam.source, beam_transfer.src_coords_cm);
        ARR3VECT(beam.direction, beam_transfer.direction);
        ARR3VECT(beam.isocenter, beam_transfer.iso_coords_cm);
        ARR2VECT(beam.fmap_size, beam_transfer.fmap_dims);
        ARR2VECT(beam.beamlet_size, beam_transfer.beamlet_size);
    }
    beam.sad = length(beam.isocenter - beam.source);

    // read standalone attributes
    {
        auto att = h5group.openAttribute("isocenter_type");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        beam.set_isocenter_type(buf);
    }
    {
        auto att = h5group.openAttribute("isocenter_location");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        beam.set_isocenter_location(buf);
    }
    {
        auto att = h5group.openAttribute("orientation_type");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        beam.set_orientation_type(buf);
    }
    {
        auto att = h5group.openAttribute("fmap_weights");
        hsize_t N = att.getSpace().getSimpleExtentNpoints();
        std::unique_ptr<float[]> temp(new float[N]);
        att.read(H5::PredType::NATIVE_FLOAT, temp.get());
        beam.fluence_map = std::vector<float>(&temp[0], &temp[N]);
    }
    return true;
}
int BEAM::_writeToHDF5(H5::Group& h5group) const { // Write beam data to hdf5 group as attributes
    H5::DataSpace scalarspace;

    { // write isocenter_type as string attribute
        std::string isocenter_type = get_isocenter_type();
        H5::StrType str_t(H5::PredType::C_S1, isocenter_type.size()+1);
        auto att = h5group.createAttribute("isocenter_type", str_t, scalarspace);
        att.write(str_t, isocenter_type);
    }
    { // write isocenter_location as string attribute
        std::string isocenter_location = get_isocenter_location();
        H5::StrType str_t(H5::PredType::C_S1, isocenter_location.size()+1);
        auto att = h5group.createAttribute("isocenter_location", str_t, scalarspace);
        att.write(str_t, isocenter_location);
    }
    { // write orientation_type as string attribute
        std::string orientation_type = get_orientation_type();
        H5::StrType str_t(H5::PredType::C_S1, orientation_type.size()+1);
        auto att = h5group.createAttribute("orientation_type", str_t, scalarspace);
        att.write(str_t, orientation_type);
    }
    { // write fluence map weights as attribute
        hsize_t fmap_dims[] = { fmap_size.x, fmap_size.y };
        H5::DataSpace fmap_ds(2, fmap_dims);

        auto att = h5group.createAttribute("fmap_weights", H5::PredType::IEEE_F32LE, fmap_ds);
        att.write(H5::PredType::NATIVE_FLOAT, fluence_map.data());
    }

    // write complex BEAM structure
    {
        COMP_BEAM_T beam_transfer;
        beam_transfer.uid = uid;
        beam_transfer.gantry_rot_rad = azimuth;
        beam_transfer.couch_rot_rad  = zenith;
        beam_transfer.coll_rot_rad   = coll;
        VECT3ARR(beam_transfer.src_coords_cm, source);
        VECT3ARR(beam_transfer.direction, direction);
        VECT3ARR(beam_transfer.iso_coords_cm, isocenter);
        VECT2ARR(beam_transfer.fmap_dims, fmap_size);
        VECT2ARR(beam_transfer.beamlet_size, beamlet_size);

        // get compound representations
        H5::CompType beam_mem_t = BEAM::getMemCompoundType();
        H5::CompType beam_file_t = BEAM::getFileCompoundType();

        auto attr = h5group.createAttribute("beam_specs", beam_file_t, scalarspace);
        attr.write(beam_mem_t, &beam_transfer);
    }
    return true;
}

void BEAM::reconfigure() {
    // Recalculate dependent parameters (source, direction, angles, ...) based on current spec
    switch (orientation_type) {
        case BEAM::ORIENT_T::ANGLE:
            // resolve ambiguities
            azimuth = fmodf(azimuth,(2*PI));
            zenith =  fmodf(zenith, (2*PI));
            // Use simplest setup possible (could affect interpretation of collimator angles)
            if (!closeto(zenith, 0.f) && (closeto(azimuth, 0.f) || closeto(azimuth, PI))) {
                zenith = 0.f;
                std::cout << set_color(COLOR::YELLOW) << "WARNING: Beam zenith (couch) angle was set to 0 deg to resolve ambiguity with azimuth (gantry) angle of "<<
                    deg_azimuth()<<" deg."<<set_color() << std::endl;
            }
            // // enforce couch range: [-90, 90]
            // if (fabsf(zenith) > PI/2.f) {
            //     zenith = -fmodf(zenith, PI/2.f);
            //     azimuth = -azimuth;
            // }

            source = BEAM::calc_source_from_angles(azimuth, zenith, isocenter, sad);
            direction = BEAM::calc_dir_from_source(isocenter, source);
            break;

        case BEAM::ORIENT_T::DIRECTION:
            // Calculate beam "angles" from this direction
            float2 angles = BEAM::calc_angles_from_dir(direction, isocenter, sad);
            azimuth = angles.x;
            zenith  = angles.y;
            source = BEAM::calc_source_from_dir(direction, isocenter, sad);
            break;
    }
}

// OP OVERLOADS
std::ostream& operator<<(std::ostream& os, const BEAM& obj) {
    // stdout print info for beam object
    std::string iso_type = obj.get_isocenter_type();
    std::string orient_type = obj.get_orientation_type();

    char outstring[200];
    sprintf(outstring,
            "[gantry=%6.1f; couch=%5.1f; coll=%5.1f || sad=%5.1f; src=(%6.1f, %6.1f, %6.1f); dir=[%4s](%6.1f, %6.1f, %6.1f); iso=[%6s](%6.1f, %6.1f, %6.1f)]",
            obj.deg_azimuth(),
            obj.deg_zenith(),
            obj.deg_coll(),
            obj.sad,
            obj.source.x, obj.source.y, obj.source.z,
            orient_type.c_str(), obj.direction.x, obj.direction.y, obj.direction.z,
            iso_type.c_str(), obj.isocenter.x, obj.isocenter.y, obj.isocenter.z
           );
    os << std::string(outstring);
    return os;
}

int get_line_count(const std::string& fname) {
    int count = 0;
    std::string line;

    std::ifstream file(fname);
    while (getline(file, line)) { count++; }

    return count;
}

// # of fields in each spec module (including identifier)
#define BEAM_SPEC_MANDATORY 2
#define BEAM_SPEC_OPT_ISO   4
#define BEAM_SPEC_OPT_DIR   4
#define BEAM_SPEC_OPT_SAD   2
int load_beam_list( std::vector<BEAM>& beams, std::string filepath, int requested_beam_count, int verbose ) {
    /* BEAM spec format:
     * Each line should contain the mandatory fields "azimuth_deg zenith_deg" as "%f %f %f"
     *   additionally, optional fields can also be provided as described in the following lines:
     *
     *   OPTIONAL 1: manual isocenter spec - "iso: %f %f %f"
     *   OPTIONAL 2: manual sad spec       - "sad: %f"
     *   OPTIONAL 3: manual beam dir spec  - "dir: %f %f %f"
     */

    // if requested_beam_count == -1, infer #beams by searching until empty line is reached
    if ( requested_beam_count <= 0 ) { requested_beam_count = get_line_count(filepath); }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Cannot open beamlist: \""<<filepath<<"\" for reading" << std::endl;
        std::cout << "failed with error ("<<errno<<"): " << std::strerror(errno) << std::endl;
        return false;
    }

    beams.resize(requested_beam_count, BEAM{});
    for (int b=0, currentline=1; b<requested_beam_count; ++b, ++currentline) {
        beams[b].sad = DEFAULT_SAD;

        if (verbose >= 2) { std::cout << "Parsing spec for beam #"<<b+1<<":"<<std::endl; }

        std::string line;
        std::getline(file, line);
        if (file.eof()){
            beams.resize(b);
            if (verbose>=2) { std::cout << "breaking early - end of file reached"<<std::endl; }
            break;
        }
        if (line.empty()) {
            b--;
            if (verbose>=2) { std::cout << "ignoring empty line: "<<currentline<<std::endl; }
            continue;
        }
        // check for comment line, skip
        if (is_comment_string(line, '#')) {
            b--;
            if (verbose>=2) { std::cout << "ignoring comment on line: "<<currentline<<std::endl; }
            continue;
        }

        // split line into fields
        std::vector<std::string> fields {};
        tokenize_string(line, fields, " ");
        int N = fields.size(); // running count of unprocesses fields
        int ptr = 0;  // position in fields vect for next unprocessed field

        // Check for mandatory fields
        if ((N-=BEAM_SPEC_MANDATORY) < 0) {
            std::cout << "ERROR: invalid beam spec detected on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
            std::cout << "  must conform to \"%f %f\"" << std::endl;
            return false;
        }
        // assign mandatories
        beams[b].azimuth = (std::stof(fields[ptr++]) * PI)/180.0f;
        beams[b].zenith  = (std::stof(fields[ptr++]) * PI)/180.0f;
        beams[b].orientation_type = BEAM::ORIENT_T::ANGLE;

        // check optionals
        while (N>0 && ptr<fields.size()) {
            if (verbose>=2) {
                std::cout << "  --> remaining fields: {";
                for (int iii=ptr; iii<N+ptr; ++iii) {
                    if (iii!=ptr) { std::cout << ", "; }
                    std::cout << "\""<<fields[iii]<<"\""; }
                std::cout << "}" << std::endl;
            }
            if (fields[ptr] == "iso:") {
                ptr++;
                if ((N-=BEAM_SPEC_OPT_ISO) < 0) {
                    std::cout << "ERROR: invalid beam spec for OPTIONAL-ISOCENTER on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
                    std::cout << "  must conform to \"iso: %f %f %f\"" << std::endl;
                    return false;
                }
                // isocenter specified in file in [units: cm]
                beams[b].isocenter.x = std::stof(fields[ptr++]);
                beams[b].isocenter.y = std::stof(fields[ptr++]);
                beams[b].isocenter.z = std::stof(fields[ptr++]);
                beams[b].isocenter_type = BEAM::ISO_T::MANUAL;
                if (verbose >= 2) { printf("  found iso=(%f, %f, %f)\n", beams[b].isocenter.x,beams[b].isocenter.y,beams[b].isocenter.z); }
                continue;
            }

            else if (fields[ptr] == "sad:") {
                ptr++;
                if ((N-=BEAM_SPEC_OPT_SAD) < 0) {
                    std::cout << "ERROR: invalid beam spec for OPTIONAL-SAD on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
                    std::cout << "  must conform to \"sad: %f\"" << std::endl;
                    return false;
                }
                // sad specfied in file in [units: cm]
                beams[b].sad = std::stof(fields[ptr++]);
                if (verbose >= 2) { printf("  found sad=%f\n", beams[b].sad); }
                continue;
            }

            else if (fields[ptr] == "dir:") {
                ptr++;
                if ((N-=BEAM_SPEC_OPT_DIR) < 0) {
                    std::cout << "ERROR: invalid beam spec for OPTIONAL-DIRECTION on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
                    std::cout << "  must conform to \"dir: %f %f %f\"" << std::endl;
                    return false;
                }
                // beam direction specified in file as unitless vector (will be normalized automatically)
                float3 dir;
                dir.x = std::stof(fields[ptr++]);
                dir.y = std::stof(fields[ptr++]);
                dir.z = std::stof(fields[ptr++]);
                beams[b].direction = normalize(dir);
                beams[b].orientation_type = BEAM::ORIENT_T::DIRECTION;

                if (verbose >=2) { printf("  found dir=(%f, %f, %f)\n", beams[b].direction.x, beams[b].direction.y, beams[b].direction.z); }
                continue;
            }

            // test for collimator rotation
            else if (dcio::is_number(fields[ptr])) {
                if ((N-=1) < 0) {
                    std::cout << "ERROR: invalid beam spec detected on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
                    std::cout << "  must conform to \"%f %f %f\"" << std::endl;
                    return false;
                }
                float coll_deg = std::stof(fields[ptr++]);
                beams[b].coll = (coll_deg * PI)/180.f;
                if (verbose >= 2) { printf("  found coll=%f\n", coll_deg ); }
                continue;
            }


            // try again starting at next field
            ++ptr;
            --N;
        }
        if (verbose && ptr < (int)fields.size()) {
            std::cout << "WARN: ignoring fields "<<ptr+1<<"-->end on line "<<currentline<<" of \""<<filepath<<"\"" << std::endl;
        }

        // Finish up
        beams[b].reconfigure();
        if (verbose>=2) { std::cout << "  Beam..."<<std::setw(4)<<std::setfill('.')<<b+1<<": "<<beams[b] << std::endl; }
    }
    return true;
}
int load_omni_beam_list( BEAM* beams, int beam_count, int verbose )
{
    std::ostringstream beamfile;
    beamfile << Paths::Instance()->temp_dir() << "/omni_beam_list.txt";
    FILE *beam_in;
    if ( (beam_in = fopen(beamfile.str().c_str(),"r")) == NULL ) {
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

        beams[b].fluence_map = std::vector<float>(beams[b].fmap_size.x * beams[b].fmap_size.y);
        if (verbose) { printf("  "); }
        load_fluence_map( beams[b].fluence_map, beams[b].fmap_size, beams[b].fmap_fname );
        if (verbose) { printf("\n"); }
    }
    fclose(beam_in);
    return 0;
}
int write_omni_beam_list( std::vector<BEAM>& beams, int beam_count, bool verbose )
{
    std::ostringstream beamfile;
    beamfile << Paths::Instance()->temp_dir() << "/omni_beam_list.txt";
    FILE *beam_out;
    if ( (beam_out = fopen(beamfile.str().c_str(),"w")) == NULL )
    {
        printf("\"write_omni_beam_list()\" failed with error (%d): %s\n", errno, std::strerror(errno));
        printf("Cannot open beam list for writing!\n");
        return -1;
    }
    fprintf(beam_out,"%d\n",beam_count);
    for (int b=0; b<beam_count; b++)
    {
        std::string isocenter_type;
        fprintf(beam_out,"%d %f %f %f %f %s %s %f %f %f %f %f %f %s %f %f %f %f %f %d %d %s\n",
                        b, beams[b].deg_azimuth(), beams[b].deg_zenith(), beams[b].deg_coll(), beams[b].sad,
                        beams[b].get_isocenter_type().c_str(),
                        beams[b].get_isocenter_location().c_str(),
                        beams[b].isocenter.x, beams[b].isocenter.y, beams[b].isocenter.z,
                        beams[b].source.x, beams[b].source.y, beams[b].source.z,
                        beams[b].get_orientation_type().c_str(),
                        beams[b].direction.x, beams[b].direction.y, beams[b].direction.z,
                        beams[b].beamlet_size.x, beams[b].beamlet_size.y,
                        beams[b].fmap_size.x, beams[b].fmap_size.y,
                        beams[b].fmap_fname.c_str() );
    }
    fprintf(beam_out,"\n#beam_id  azimuth zenith coll sad | iso_type iso_loc isocenter.x isocenter.y isocenter.z | source.x source.y source.z | orient_type direction.x direction.y direction.z | beamlet_size.x beamlet_size.y fluence_dim.x fluence_dim.y | fmap_fname\n");
    fclose(beam_out);
    return 0;
}
