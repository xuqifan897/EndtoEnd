#ifndef __BEAM_H__
#define __BEAM_H__

#include <iostream>
#include <string>
#include <vector>

#include "helper_cuda.h"
#include "helper_math.h"
#include "math_constants.h"
#include "./macros.h"

#include "H5Cpp.h"

// beam description
class BEAM
{
    public:
        enum class ORIENT_T {
            ANGLE,
            DIRECTION,
            UNSPEC
        };
        enum class ISO_T {
            UNSPEC,
            PTV_CENTROID,
            MANUAL
        };
        enum class ISO_LOC_T {
            IN_PTV,
            OUT_PTV,
            UNSPEC
        };

        uint   uid;             // unique index that maps beam back to line in beam_list.txt
        float  azimuth;         // Gantry Angle     [rad]
        float  zenith;          // Couch Angle      [rad]
        float  coll;            // Collimator Angle [rad]
        float3 source;          // x-ray src coords [cm]
        float3 direction;       // norm. direction vector
        ORIENT_T orientation_type;   // Method used to select orientation
        float3 isocenter;       // isocenter coords [cm]
        ISO_T  isocenter_type;  // Method used to select isocenter
        ISO_LOC_T isocenter_location; // describe in or out of PTV
        float  sad;             // source-axis-distance [cm]

        uint2  fmap_size;       // fluence map dims
        float2 beamlet_size;    // beamlet size [cm]
        std::string fmap_fname; // optimization results / binary mask
        std::vector<float> fluence_map; // linear-storage beamlet weights

        CUDEV_FXN float deg_azimuth() const { return 180.0f * azimuth / CUDART_PI_F; }
        CUDEV_FXN float deg_zenith()  const { return 180.0f * zenith / CUDART_PI_F; }
        CUDEV_FXN float deg_coll()    const { return 180.0f * coll / CUDART_PI_F; }
        static CUDEV_FXN float3 calc_source_from_angles(float gantry_rot_rad, float couch_rot_rad, float3 iso, float sad);
        static CUDEV_FXN float2 calc_angles_from_dir(float3 dir, float3 iso, float sad);
        static CUDEV_FXN float3 calc_source_from_dir(float3 dir, float3 iso, float sad);
        static CUDEV_FXN float3 calc_dir_from_source(float3 iso, float3 source);

        CUDEV_FXN float3 calc_source_from_angles() const;
        CUDEV_FXN float2 calc_angles_from_dir() const;
        CUDEV_FXN float3 calc_source_from_dir() const;
        CUDEV_FXN float3 calc_dir_from_source() const;

        void set_isocenter_type(const std::string& iso_type);
        std::string get_isocenter_type() const;
        void set_isocenter_location(const std::string& iso_loc);
        std::string get_isocenter_location() const;
        void set_orientation_type(const std::string& orient_type);
        std::string get_orientation_type() const;
        void reconfigure();

        static H5::CompType getFileCompoundType();
        static H5::CompType getMemCompoundType();
        static int _readFromHDF5(BEAM& beam, H5::Group& h5group);
        int _writeToHDF5(H5::Group& h5group) const;

        friend std::ostream& operator<<(std::ostream& os, const BEAM& obj);
    private:
        // POD describing compound HDF5 datatype used in storing beam metadata
        struct COMP_BEAM_T {
            ushort      uid;              // Unique ID (0-based)

            float       gantry_rot_rad;   // Gantry Angle     [rad]
            float       couch_rot_rad;    // Couch Angle      [rad]
            float       coll_rot_rad;     // Collimator Angle [rad]
            float       src_coords_cm[3]; // x-ray src coords [cm]
            float       direction[3];     // norm. direction vector
            float       iso_coords_cm[3]; // isocenter coords [cm]

            uint        fmap_dims[2];     // fluence map dims
            float       beamlet_size[2];  // beamlet size [cm]
        };

};
std::ostream& operator<<(std::ostream& os, const BEAM& obj);

int load_beam_list( std::vector<BEAM>& beams, std::string filepath, int requested_beam_count=-1, int verbose=false);
int load_omni_beam_list( BEAM* beams, int beam_count, int verbose=0 );
int write_omni_beam_list( std::vector<BEAM>& beams, int beam_count, bool verbose=false );

#endif // __BEAM_H__
