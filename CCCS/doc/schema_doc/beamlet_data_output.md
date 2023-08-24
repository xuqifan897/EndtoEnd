## Data Hierarchy:
Description of HDF5 group/attribute/dataset hierarchy with valid field names, datatypes, and element counts (`!|-` indicates a conditional field only added when ROI-masking/row-reduction is performed)
```c++
GROUP "/" {
|-  GROUP "filetype" {
|   |-  ATTRIBUTE "ftmagic"               { [1] H5T_STD_U8LE = 43 }
|   |-  ATTRIBUTE "ftversionmajor"        { [1] H5T_STD_U8LE = 1 }
|   |-  ATTRIBUTE "ftversionminor"        { [1] H5T_STD_U8LE = 2 }
|   }
| 
|+  GROUP "calc_specs" {
|   |-  ATTRIBUTE "N_beams"               { [1] H5T_STD_U16LE }
|   |-  ATTRIBUTE "calc_bbox_size"        { H5T_ARRAY { [3] H5T_STD_U16LE } }
|   |-  ATTRIBUTE "calc_bbox_start"       { H5T_ARRAY { [3] H5T_STD_U16LE } }
|   |-  ATTRIBUTE "penumbra_cm"            { [1] H5T_IEEE_F32LE }
|   |-  ATTRIBUTE "convlat_cm"            { [1] H5T_IEEE_F32LE }
|   |-  ATTRIBUTE "convstep_cm"           { [1] H5T_IEEE_F32LE }
|   |-  ATTRIBUTE "dicom_start_cm"        { H5T_ARRAY { [3] H5T_IEEE_F32LE } }
|   |-  ATTRIBUTE "full_dicom_size"       { H5T_ARRAY { [3] H5T_STD_U16LE } }
|   |-  ATTRIBUTE "nphi"                  { [1] H5T_STD_U16LE }
|   |-  ATTRIBUTE "nradii"                { [1] H5T_STD_U16LE }
|   |-  ATTRIBUTE "ntheta"                { [1] H5T_STD_U16LE }
|   |-  ATTRIBUTE "sparsity_thresh"       { [1] H5T_IEEE_F32LE }
|   |-  ATTRIBUTE "voxel_size_cm"         { H5T_ARRAY { [3] H5T_IEEE_F32LE } }
|   |-  ATTRIBUTE "beam_spectrum"         { [str_len] H5T_STRING }
|   |-  ATTRIBUTE "target_structure"      { [str_len] H5T_STRING }
|  !|-  ATTRIBUTE "roi_order"             { [N_roi [variable] ] H5T_STRING }
|  !|-  ATTRIBUTE "row_block_capacities"  { [N_roi] H5T_STD_U64LE }
|   }
|
|+  GROUP "beams" {
|   |+  GROUP "data" {
|   |   |+  GROUP "beam_00000" {
|   |   |   |+  GROUP "beamlet_00000" {
|   |   |   |   |-  ATTRIBUTE "N_coeffs"              { [1] H5T_STD_U32LE }
|   |   |   |   |-  ATTRIBUTE "beamlet_uid"           { [3] H5T_STD_U16LE }
|   |   |   |   |-  ATTRIBUTE "perc_dropped"          { [1] H5T_IEEE_F32LE DATASPACE }
|   |   |   |   |-  ATTRIBUTE "perc_nonzero"          { [1] H5T_IEEE_F32LE DATASPACE }
|   |   |  !|   |-  ATTRIBUTE "row_block_sizes        { [N_roi] H5T_STD_U64LE } 
|   |   |   |   |-  DATASET "coeffs"                  { [N_coeffs] H5T_IEEE_F32BE }
|   |   |   |   |-  DATASET "lindex"                  { [N_coeffs] H5T_STD_U32LE }
|   |   |   |   }
|   |   |   |   
|   |   |   :+  
|   |   |   }  
|   |   |   
|   |   :+  
|   |   }  
|   |   
|   |+  GROUP "metadata" {
|   |   |+  GROUP "beam_00000" {
|   |   |   |-  ATTRIBUTE "N_beamlets"            { [1] H5T_STD_U16LE }
|   |   |   |-  ATTRIBUTE "beam_uid"              { [1] H5T_STD_U16LE }
|   |   |   |-  ATTRIBUTE "isocenter_type"        { [str_len] H5T_STRING }
|   |   |   |-  ATTRIBUTE "beam_specs"            { H5T_COMPOUND {
|   |   |   |         [1] H5T_STD_U16LE  "uid";
|   |   |   |         [1] H5T_IEEE_F32LE "gantry_rot_rad";
|   |   |   |         [1] H5T_IEEE_F32LE "couch_rot_rad";
|   |   |   |         [1] H5T_IEEE_F32LE "coll_rot_rad";
|   |   |   |         H5T_ARRAY { [3] H5T_IEEE_F32LE } "src_coords_cm";
|   |   |   |         H5T_ARRAY { [3] H5T_IEEE_F32LE } "direction";
|   |   |   |         H5T_ARRAY { [3] H5T_IEEE_F32LE } "iso_coords_cm";
|   |   |   |         H5T_ARRAY { [2] H5T_STD_U32LE }  "fmap_dims";
|   |   |   |         H5T_ARRAY { [2] H5T_IEEE_F32LE } "beamlet_size_cm";
|   |   |   |     }}
|   |   |   |-  ATTRIBUTE "fmap_weights"          { [beam_specs.fmap_dims] H5T_IEEE_F32LE }
|   |   |   |   
|   |   |   :+ 
|   |   |   }
|   |   |
|   |   :+
|   |   }
|
|
}
```

## Data Definitions
All data structures are documented in C-struct format. Explicit storage types, and bit-order (endianness) are maintained in HDF5 files and any capable HDF5 library should handle data interpretation and type conversions automatically.

### Patient/Problem Attributes [`/calc_specs/`]:
```c++
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
    float       delr_cm;              // CCCS convolution stepsize in [cm]
    float       kernel_extent_cm;     // maximum radial distance used in CCCS [cm]
    std::string beam_spectrum;        // name of beam spectrum file (describes beam energy)
    std::string target_structure;     // name of structure to which dose is delivered

    // Added after reduction of column
    std::vector<std::string> roi_order{};            // ordered list of roi_names for Row-Blocks in reduced A-matrix
    std::vector<uint64_t>    row_block_capacities{}; // ord. list of row block capacities matching roi_order
};
```
### Beam Attributes [`/beam_###/`]:
```c++
struct HEADER_BEAM {
    ushort      beam_uid;         // UID for beam (index into patient_header.beam_list)
    BEAM        beam_specs;       // parameters specifying beam (angles, source location, ...)
    ushort      N_beamlets;       // number of active beamlets in this beam
};
struct BEAM
{
    uint16_t uid               // UID for beam (index into patient_header.beam_list)
    float    gantry_rot_rad;   // Gantry Angle     [rad]
    float    couch_rot_rad;    // Couch Angle      [rad]
    float    coll_rot_rad;     // Collimator Angle [rad]
    float3   src_coords_cm;    // x-ray src coords [cm]
    float3   direction;        // norm. direction vector
    float3   iso_coords_cm;    // isocenter coords [cm]

    uint2    fmap_dims;        // fluence map dims
    float2   beamlet_size_cm;  // beamlet size [cm]
};
enum class ISO_T {
    UNSPEC,                    // no type specified.  as_str: "unspec"
    PTV_CENTROID,              // auto from ptv.      as_str: "ptv"
    MANUAL                     // spec. in beamlist.  as_str: "man"
};

```
### Beamlet Attributes [`/beam_###/beamlet_###/`]:
```c++
struct HEADER_BEAMLET {
    ushort      beamlet_uid;      // UID for beamlet (linearized index into 2D fluence map)
    uint64_t    N_coeffs;         // number of non-zero coeffs in column

    // added after reduction of column
    std::vector<uint64_t> row_block_sizes;                      // list indicating true non-zeros per row-block after reduction
};
```
### Sparse Dataset: [`/beam_###/beamlet_###/{lindex,coeffs}`]:
On disk, separate arrays are stored side-by-side (see [above](#data-hierarchy)), in memory these are collected into a C-structure
```c++
struct KeyValPairs {
    std::vector<uint64_t>  keys;  // linearized index into full dicom volume
    std::vector<float>     vals;  // dose coefficient value
};
```

## Appendix A: CUDA Vector Types
In-memory format of CUDA vectors:
```c++

struct __device_builtin__ float2
{
    float x, y;
};
struct __device_builtin__ float3
{
    float x, y, z;
};
struct __device_builtin__ uint2
{
    unsigned int x, y;
};
struct __device_builtin__ uint3
{
    unsigned int x, y, z;
};
```
On-disk format (converted to N-element arrays):
```c++
[ x, y, z ]
```
