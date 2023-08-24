# Introduction:
This page describes the output file format to which the beam angle selection and fluence map optimizer (FMO) conforms. This format is essentially a subset of the beamlet-based dose calculation output format described [here](https://github.com/radiasoft/rs4pi/wiki/DoseCalc-Data-Output-Specification) in that it describes the beam specifications but excludes beamlet level information.

This file format is fully compatible with dosecalc-preprocess and can be used with the `--fmaps=` cli arg to load selected beams and optimized fluence intensities in preparation for further beamlet coefficient calculation (dosecalc-beamlet) or final dose calculation (dosecalc_beam).

## Data Hierarchy:
Description of HDF5 group/attribute/dataset hierarchy with valid field names, datatypes, and element counts:
```c++
GROUP "/" {
|-  GROUP "filetype" {
|   |-  ATTRIBUTE "ftmagic"               { [1] H5T_STD_U8LE = 44 }
|   |-  ATTRIBUTE "ftversionmajor"        { [1] H5T_STD_U8LE = 1 }
|   |-  ATTRIBUTE "ftversionminor"        { [1] H5T_STD_U8LE = 1 }
|   }
|
|+  GROUP "beams" {
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
|   |
|
}
```

## Additional Information:
Refer to the extended documentation in [DoseCalc Data Output Specification](https://github.com/radiasoft/rs4pi/wiki/DoseCalc-Data-Output-Specification)

