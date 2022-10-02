# EndtoEnd
Continuous optimization of photon dose radiotherapy

## Background
  X-ray radiation therapy aims to treat tumors or other lesions while sparing normal tissue. For this purpose, treatment planning is used to select the optimal beam orientations and the fluence map through optimization. Previous methods used linear and discrete model, in which the fluence map was subdivided into pixels or beamlets, and the dose was calculated as:
  $$d_i = \sum_jB_{i,j}w_j,$$
  where $d_i$ denotes the dose value at the $i^{th}$ voxel, $B_{i,j}$ is the dose value of the $j^{th}$ beamlet at the $i^{th}$ voxel, and $w_j$ the weight of the $j^{th}$ beamlet. Or in matrix form:
  $$\vec{d}=B\vec{w}.$$
  To achieve this, the dose matrix $B$ had to be pre-calculated and stored. Beam selection was achieved by imposing a sparsity term with a tunable weight controlling the number of beams selected.  
  However, due to the memory constraint and computation burden, the fluence map was subdivided with low granularity, and only a limited number of candidate beams were selected, which potentially resulted in suboptimal results.
  
  With the help of the massive parallelism of GPU, it is possible to compute the dose on the fly, without the need for pre-computing and storage. In this project, utilizing the linear relationship between the fluence map and the dose, the gradient to fluence map can be simpliy calculated using the chain rule, and the orientation parameters were updated using perturbation.

## Methods overview
  Following the [FCBB](FCBB_link) methods, there are three steps to compute the dose from the fluence map for a single beam:
  * Fluence map convolution
  * Beam's eye view (BEV) dose calculation
  * Patient volume coordinate system (PVCS) dose calculation

  Fluence map convolution addresses the dose spread caused by secondary particles. BEV dose calculation calculates the dose in BEV view, and PVCS dose calculation transforms the dose in BEV to PVCS. All the three steps are linear operations, which makes gradient computation feasible.

## Quick start
From here, we assume the project folder is at root. To configure and compile the code, run
```
user@host:~$ bash configure.sh
user@host:~$ bash build.sh
```
This produces five executables in ./build directory:
* `dose_calculation`: to calculate the dose of each input beams. An example execution script is given in `./examples/dose_calculation.sh`. Before running it, the user needs to change the argument list accordingly. Then run

        user@host:~$ bash ./examples/dose_calculation.sh
  To get the instruction of the parameter specification, run

        user@host:~$ ./build/dose_calculation --help

  This works for other executables as well, just replace the executable name.

* `optimize_stationary`: to optimize the fluence map, without changing the beam angles. The loss function only involves the dose consistency, without considering smoothness. To run,

        user@host:~$ bash ./examples/optimize_stationary.sh

* `optimize_stationary_smoothness`: to optimize the fluence map, without changing the beam angles. The loss function is a combination of dose consistency and fluence map smoothness. To run,

        user@host:~$ bash ./examples/optimize_stationary_smoothness.sh
* `optimize_dynamic`: to collectively optimize fluence map and beam angles. To run,

        user@host:~$ bash ./examples/optmize_dynamic.sh
* `optimize_dynamic_random.sh`: to collectively optimize fluence map and beam angles, while integrating some randomness to prevent from being trapped in local minimum (which proved to faile). To run,

        user@host:~$ bash ./examples/optimize_dynamic_random.sh

## Arguments explaination
Here we address the details of some arguments.
* `phantom-path`: it is a path to a phantom, which is stored in binary format. To visualize the phantom, the user can use numpy to load it and reshape it into the desired dimension.

        import numpy as np
        phantom = np.fromfile(${phantom-path}, dtype=np.float32)
        phantom = np.reshape(phantom, )

[FCBB_link]: https://pubmed.ncbi.nlm.nih.gov/21081826/