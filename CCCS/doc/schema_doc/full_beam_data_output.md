# Introduction:
This page describes the data file format produced by the final dose calculator (mgcs_omnidose_full).

## Data Hierarchy:
```c++
GROUP "/" {
|+ GROUP "filetype" {
|  |-  ATTRIBUTE "ftmagic"               { [1] H5T_STD_U8LE = 42 }
|  |-  ATTRIBUTE "ftversionmajor"        { [1] H5T_STD_U8LE = 1 }
|  |-  ATTRIBUTE "ftversionminor"        { [1] H5T_STD_U8LE = 0 }
|  }
| 
|+ DATASET "dose" {
|  |  DATATYPE  H5T_IEEE_F32LE
|  |  DATASPACE  SIMPLE { ( 256, 256, 256 ) / ( 256, 256, 256 ) }
|  |
|  |- ATTRIBUTE "dicom_start_cm" { [3] H5T_IEEE_F32LE }
|  |- ATTRIBUTE "voxel_size_cm" { [3] H5T_IEEE_F32LE }
|  }
}
```
 Note: the data is in ZYX order (x is fastest index, z is slowest index). Dataset size reflects ZYX ordering
