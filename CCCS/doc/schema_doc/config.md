Config.json Schema
==================
|                        |             |
|------------------------|-------------|
| **Schema Version**     | 1.3         |
| **Structure**          | text/json   |
| **Last Modified By**   | Ryan Neph   |
| **Last Modified Date** | Apr 29 2020 |


------------------------------------------------------------------------

| **option**    | **acceptable data types** | **Description**                                                                                                |
|---------------|---------------------------|----------------------------------------------------------------------------------------------------------------|
| voxsize       | float                     | isotropic voxelsize assumed for production of dose volume (units: cm)                                          |
| convlat       | float                     | Lateral ray spacing used during CCCS convolution (units: cm)                                                   |
| convstep      | float                     | Longitudinal step size used during CCCS convolution (units: cm)                                                |
| beamlet-size  | float; [float, float]     | 2D size of beamlet at isocenter (units: cm)                                                                    |
| fmap-dims     | float; [float, float]     | 2D dimensions of fluence map (MLC)                                                                             |
| penum         | float                     | expansion of dose calculation beyond terma bounding box (units: cm) [omnidose-full only]                       |
| kernel-extent | float                     | approximation parameter terminating CCCS convolution after dose spread at this radial distance (units: cm)     |
| ntheta        | uint                      | rebinning of dose kernel into fewer dose spread angles                                                         |
| nphi          | uint                      | kernel axial rotation count                                                                                    |
| max-rev-size  | [uint, uint, uint]        | size of static convolution space on GPU (affects overall GPU RAM usage)                                        |
| verbose       | bool                      | enable verbose output                                                                                          |
| timing        | bool                      | enable timing output                                                                                           |
| nbeams        | uint                      | limit number of beams to read from _beamlist_                                                                  |
| spec          | string                    | specify name of beam spectrum definition file (w/o path)                                                       |
| target        | string                    | use first structure containing this substring as target                                                        |
| target-exact  | string                    | use first structure exactly matching this string as target                                                     |
| bbox-roi      | string                    | use first structure containing this substring as dose calculation volume                                       |
| beamlist      | string                    | path to file containing specification of beams to calculate                                                    |
| fmaps         | string                    | path to file containing beam orientation specifications and fluence map intensities                            |
| ctlut         | string                    | path to file containing CT# to material mass density table.                                                    |
