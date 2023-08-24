#ifndef __IO_DATA_STRUCT_H__
#define __IO_DATA_STRUCT_H__

#include <string>
#include <vector>
#include <helper_cuda.h>
#include <helper_math.h>

#include "./beam.h"

// POD for sparsified coefficients
struct KeyValPairs {
    std::vector<uint64_t>  keys;  // linearized index into full dicom volume
    std::vector<float>     vals;  // dose coefficient value
    size_t size() const { return keys.size(); }
};
// Encapsulating storage describing a full column of the A-matrix, supporting [non-]reduced structures
struct SparseData {
    KeyValPairs kvp;
    size_t size() const { return kvp.size(); }
    float perc_nonzero; // percentage of calc_bbox that remains as nonzero (above thresh)
    float perc_dropped; // percentage of calc_bbox that has been zero'd (below thresh)

    // extra data added after row reduction/reordering
    std::vector<uint64_t> row_block_sizes;
    bool reduced() { return !row_block_sizes.empty(); }
    float perc_reducedrop; // percentage of calc_bbox that has been excluded due to roi masking
};
// POD describing storage of meta to HDF5 group for each beamlet (column)
struct HEADER_BEAMLET {
    ushort      beamlet_uid;      // UID for beamlet (linearized index into 2D fluence map)
    uint64_t    N_coeffs;         // number of non-zero coeffs in column

    // added after reduction of column
    std::vector<uint64_t> row_block_sizes;                      // list indicating true non-zeros per row-block after reduction
    bool  reduced() const { return !row_block_sizes.empty(); }; // flag indicating whether reduction occured
};
// POD describing storage of meta to HDF5 group for each beam (collection of beamlets)
struct HEADER_BEAM {
    ushort      beam_uid;         // UID for beam (index into patient_header.beam_list)
    BEAM        beam_specs;       // parameters specifying beam (angles, source location, ...)
    ushort      N_beamlets;       // number of active beamlets in this beam
};
// POD describing storage of meta to HDF5 group for each full calculation/execution
// contents are fields that only need to be stored in one instance for each output file
struct HEADER_PATIENT {
    ushort      N_beams;              // total number of beams stored in file
    float3      dicom_start_cm;       // coords of dicom start in [cm]
    uint3       full_dicom_size;      // size of full dicom volume
    float3      voxel_size_cm;        // voxel size in [cm]
    float       rev_latspacing_cm;    // Conv ray lateral spacing [cm]
    float       rev_longspacing_cm;   // Conv longitudinal step size [cm]
    uint3       bbox_start;           // indices of start of dose calc_bbox
    uint3       bbox_size;            // size of dose calc_bbox
    float       penumbra_cm;          // beam margin for dose calculation [cm]
    float       sparsity_thresh;      // Setting for sparse thresholding
    ushort      nphi;                 // number of zenithal kernel rays
    ushort      ntheta;               // number of azimuthal kernel rays
    ushort      nradii;               // number of radial samples along each kernel ray
    float       kernel_extent_cm;     // maximum radial distance used in CCCS [cm]
    std::string beam_spectrum;        // name of beam spectrum file (describes beam energy)
    std::string target_structure;     // name of structure to which dose is delivered

    // Added after reduction of column
    std::vector<std::string> roi_order{};            // ordered list of roi_names for Row-Blocks in reduced A-matrix
    std::vector<uint64_t>    row_block_capacities{}; // ord. list of row block capacities matching roi_order
    bool reduced() { return !roi_order.empty(); }  // flag indicating setting of remaining members
    bool reduced() const { return !roi_order.empty(); }  // flag indicating setting of remaining members
};

// Describes sub-volume for iterating or calculating
struct ArrayProps {
    uint3 size;                   // size of full dicom volume
    uint3 crop_size;              // size of calc_bbox
    uint3 crop_start;             // start indices of calc_bbox relative to original dicom volume size
    // FrameOfReference frame {};       // NOT IN USE YET - MIGRATE size to frame.size
    CUDEV_FXN uint3 crop_end() { return crop_start + crop_size; }
    uint nvoxels() { return crop_size.x * crop_size.y * crop_size.z; }

    static int _readFromHDF5(ArrayProps& props, H5::Group& h5group) {
        // tuple dataspaces/datatypes
        hsize_t tuple3_dims[] = { 3 };
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_UINT, 1, tuple3_dims);

        // read attributes
        {
            uint temp[3];
            auto att = h5group.openAttribute("size");
            att.read(tuple3_native_t, temp);
            ARR3VECT(props.size, temp);
        }
        {
            uint temp[3];
            auto att = h5group.openAttribute("crop_size");
            att.read(tuple3_native_t, temp);
            ARR3VECT(props.crop_size, temp);
        }
        {
            uint temp[3];
            auto att = h5group.openAttribute("crop_start");
            att.read(tuple3_native_t, temp);
            ARR3VECT(props.crop_start, temp);
        }
        return true;
    }
    int _writeToHDF5(H5::Group& h5group) const {
        // tuple dataspaces/datatypes
        H5::DataSpace scalarspace {};
        hsize_t tuple3_dims[] = { 3 };
        H5::DataSpace tuple3(1, tuple3_dims);
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_UINT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::STD_U16LE, 1, tuple3_dims);

        // write attributes
        {
            uint temp[3];
            VECT3ARR(temp, size)
            auto att = h5group.createAttribute("size", tuple3_t, scalarspace);
            att.write(tuple3_native_t, temp);
        }
        {
            uint temp[3];
            VECT3ARR(temp, crop_size)
            auto att = h5group.createAttribute("crop_size", tuple3_t, scalarspace);
            att.write(tuple3_native_t, temp);
        }
        {
            uint temp[3];
            VECT3ARR(temp, crop_start)
            auto att = h5group.createAttribute("crop_start", tuple3_t, scalarspace);
            att.write(tuple3_native_t, temp);
        }
        return true;
    }
};

// defines layout of pillars in beamlet-based dose calculation
//  all array data structures share index, beamletIdx[] maps this index to beamlet fmap linear index
struct PILLAR_GRID{
    int numBeamlets;             // number of active pillars in this fluence map
    int2 numPillars;             // number of pillars along x and z axis in pillarMap.
    int3 gridDims;               // precalculated numPillars * numVoxelPerPillar (total voxel dims of packed array)
    int3 pillarDims;             // dimension of pillars along x,y,z directions in voxels
    float3  gridSize;            // dimensions of grid in common coordinate unit
    float3  pillarSize;          // dimensions of pillars in common coordinate unit
    float2  max_beamlet_size;    // maximum size of beamlet due to divergence/magnification
    int wallThickness = 2;       // number of REV voxels composing separating walls; should be >=2 for best results
    int pillarBuffer = 2;        // number of REV voxels buffering usable pillar space from walls (garbage ends up here from interpolation near walls); Should be >=2 for best results

    // arrays
    int* beamletIdx;             // indices of actual beamlets
    float3* pillarStartCoords;   // coords of pillar box start
    float2* beamletAngles;       // sum of beam angle and beamlet divergence angle (.x: azimuth; .y: zenith)
    float3* beamletIsocenters;   // RCS coords of each beamlet isocenter (intersect of beamlet central axis and fluence map plane)
    float3* beamletRotationAnchors; // RCS-aligned offset from pillarStartCoords that defines the beamletRotation anchor point

    // convenience functions
    CUDEV_FXN inline int pillar_nvoxels() const { return pillarDims.x * pillarDims.y * pillarDims.z; }
    CUDEV_FXN inline int pillar_grid_nvoxels() const { return gridDims.x * gridDims.y * gridDims.z; }
};

// struct to hold data for size and dimension of beam's eye view fields
struct REV_DATA {
    float3 min_coords {}; // coords of first voxel in BEV
    float3 max_coords {}; // coords of last voxel in BEV
    uint3  size {};       // dimensions of this BEV volume
};

struct BLT_CONV_DATA {
    PILLAR_GRID  pgrid {};
    REV_DATA     rev {};
    unsigned int nbatches = 1;

    // cuda kernel launch params
    dim3         raytrace_grid {};
    dim3         raytrace_block {};
    unsigned int raytrace_sharedmem = 0;
};

#endif // __IO_DATA_STRUCT_H__
