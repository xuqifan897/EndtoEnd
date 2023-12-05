#ifndef __BRAIN_DEFS_H__
#define __BRAIN_DEFS_H__

#include "dosecalc_defs.h"
#include "beam.h"
#include "macros.h"

#include <vector>
#include <string>

#define _MAX_THREADS_PER_BLOCK_ 1024       // (CUDA CC 2.0+)
// #define _MAX_GRID_SIZE_         4294967295 // (CUDA CC 3.0+)
#define _MAX_GRID_SIZE_         65535 // (CUDA CC <3.0)

namespace old
{
    class CONSTANTS {
    public:
        // kernel params
        unsigned int nphi;      // kernel azimuthal angles
        unsigned int ntheta;    // kernel zenithal angles
        unsigned int nradii;    // number of radial boundaries in dose deposition kernel
        std::vector<float> conv_theta_deg;  // array of convolution theta angles indexed by "ray_index"
        std::vector<float> conv_phi_deg;    // array of convolution phi   angles indexed by "ray_index"

        // calculation params
        uint ss_factor;         // terma anti-aliasing supersample factor
        float rev_latspacing;   // convolution lateral ray spacing [unit: cm]
        float rev_longspacing;  // convolution longitudinal step size [unit: cm]
        float kernel_extent;    // dose kernel radius truncation distance [unit: cm]
        float penumbra;         // margin added to rev in fullbeam calc [unit: cm]
        float beamhard_correct; // beam hardening correction used in terma calculation

        unsigned int beam_count;
        BEAM beam;                 // 3D coordinates of beam origin, size of beam in cm
        std::string target_structure; // name of PTV
        std::string beam_spec;        // beam spectrum filename (no extension)

        // options
        bool reduce =false;               // indicator for ROI based sparse data masking and reordering (M-matrix to A-matrix inline conversion)

        // dicom properties
        float3 start;           // 3D coordinates of data origin (i.e. first voxel of the dicom volume)
        float3 voxel;           // size of voxel in mm
        uint3  size;            // dimensions of data in voxels after resampling
        CUDEV_FXN inline int nvoxels() const { return size.x*size.y*size.z; }
        CUDEV_FXN inline float3 end() const {
            return make_float3(
                    start.x + voxel.x*size.x,
                    start.y + voxel.y*size.y,
                    start.z + voxel.z*size.z
                );
        }

        // BEV properties
        uint3 max_rev_size;     // max BEV dimensions per beam orientation

        // Calculation Bounding Box Properties
        uint3 calc_bbox_start; // # of voxels from dicom origin to bounding box
        uint3 calc_bbox_size;
        CUDEV_FXN inline int bbox_nvoxels() const { return calc_bbox_size.x*calc_bbox_size.y*calc_bbox_size.z; }
        CUDEV_FXN inline float3 bbox_start_to_coords() const {
            return start + make_float3(calc_bbox_start)*voxel;
        }

        // Added for convenience and kernel precalculation
        int calc_bbox_precomp_2d_size;

        ///////////// INLINES //////////////
        // convenience functions for converting index to radian angle
        inline int get_kernel_theta_index(int ray_index) {
            return ray_index % (ntheta/2);
        }
        inline int get_kernel_phi_index(int ray_index) {
            return ray_index / (ntheta/2);
        }
        inline float get_theta_from_index(int ray_index) {
            return conv_theta_deg[ray_index]*PI/180.f;
            /* return (float)get_kernel_theta_index(ray_index, ntheta, nphi) * PI / ntheta; */
        }
        inline float get_phi_from_index(int ray_index) {
            return conv_phi_deg[ray_index]*PI/180.f;
            /* return (float)get_kernel_phi_index(ray_index, ntheta, nphi) * 2.0f * PI / nphi; */
        }
        ////////////////////////////////////

        // Flags for debugging
        bool debugREVTerma;
        bool debugREVDose;
        bool debugBEVDose;
        bool logPatientDose;
    };

    struct SHM_DATA {
        unsigned int size_data;

        std::vector<float> density;         // density data
        std::vector<float> kernel_array;    // dose deposition kernel
        std::vector<float> radial_boundary; // array of boundary radii in cm
    };
}

#endif