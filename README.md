# EndtoEnd
Direct beam angle optimization with simulated annealing based on GPU dose calculation

## Background
  X-ray radiation therapy aims to treat tumors or other lesions while sparing normal tissue. For this purpose, treatment planning is used to select the optimal beam orientations and the fluence map through optimization. Previous methods used linear and discrete model, in which the fluence map was subdivided into pixels or beamlets, and the dose was calculated as:
  $$d_i = \sum_jB_{i,j}w_j,$$
  where $d_i$ denotes the dose value at the $i^{th}$ voxel, $B_{i,j}$ is the dose value of the $j^{th}$ beamlet at the $i^{th}$ voxel when the weight of the $j^{th}$ beamlet is unity, and $w_j$ the weight of the $j^{th}$ beamlet. Or in matrix form:
  $$\vec{d}=B\vec{w}.$$
  To achieve this, the dose matrix $B$ had to be pre-calculated and stored. Beam selection was achieved by imposing a sparsity term with a tunable weight controlling the number of beams selected.  
  However, due to the memory constraint and computation burden, the fluence map was subdivided with low granularity, and only a limited number of candidate beams were selected, which potentially resulted in suboptimal results.
  
  With the help of the massive parallelism of GPU, it is possible to compute the dose on the fly, without the need for pre-computing and storage. In this project, utilizing the linear relationship between the fluence map and the dose, the gradient to fluence map can be simpliy calculated using the chain rule, and the orientation parameters are updated using simulated annealing.

## Methods overview
  Following the [FCBB] methods, there are three steps to compute the dose from the fluence map for a single beam:
  * Fluence map convolution
  * Beam's eye view (BEV) dose calculation
  * Patient volume coordinate system (PVCS) dose calculation

  Fluence map convolution addresses the dose spread caused by secondary particles. BEV dose calculation calculates the dose in BEV, and PVCS dose calculation transforms the dose in BEV to PVCS. All the three steps are linear operations, which makes gradient computation feasible. The loss function in our setting is $$L=\sum_i w_{0i}(d_i-t_{0i})^{2} + \sum_{j=1}^{OARs}(d_i-t_{ji})_{+}^{2} + \eta||Dx_k||_{1}$$ where $i$ is the index for voxels, $j$ is the index for anatomical structures ($j=0$ for PTV, $j\geq1$ for OARs), $k$ is the index for beams, $w_{ji}$ is the weight for the $i^{th}$ voxel of the $j^{th}$ anatomical structure, $d_i$ is the dose at the $i^{th}$ voxel, $t_{ji}$ is the target dose for the $i^{th}$ voxel of the $j^{th}$ anatomical structure, $(\sdot)_+$ clamps the value to non-negative, $\eta$ is the weight controlling smoothness, $D$ is the gradient operator in both axises of the fluence map, $||\sdot||_1$ takes the 1-norm. The dose and fluence map is related by $d=\sum_{k=1}^{beams}A_kx_k$.

  Intuitively, the first term penalizes the dose deviation from the target, the second term penalizes the dose greater than the target, and the third term encourages smoothness of the fluence map. It should be noted that different anatomical structures may overlap, i.e., $w_{j_1i}>0,w_{j_2i}>0$ with $j_1 \neq j_2$. In fact, if we rearrange the order of summation in the second term, we can rewrite it as a sum of terms, each of which is a quadratic term w.r.t. $d_i$ for all $i$. In other words, we can coallese $w_j$ and $t_j$ for all $j>0$ into a single weight matrix $w$ and a target matrix $t$ in practice.

## Quick start
### Data preparation
  To reproduce the results (to appear in our paper), please find the dataset through this [link](https://doi.org/10.5281/zenodo.7236751). This dataset contains the image data and anatomical annotation in dicom format. For each patient, there are two zip files. Take patient1 for example:
  * patient1_E2E\
    ├── beam_angles_annealing.txt ...... *resulting beam angles of the annealing algorithm* \
    ├── beam_angles_E2E.txt ...... *resulting beam angles of the baseline beam selection algorithm* \
    ├── beam_angles_uniform.txt ...... *equally spaced coplanar beam angles, used as initialization of the annealing algorithm* \
    ├── beam_angles_VarianIEC.csv ...... *beam_angles_E2E.txt, in VarianIEC format, i.e., gantry and couch angles* \
    ├── beam_annealing_varianIEC.csv ...... *beam_angles_annealing.txt, in VarianIEC format* \
    ├── CT.dat ...... *the binary CT data in HU (0.0 for air and 1000.0 for water), as a 32-bit float array of dimension hight\*width\*slices, in lexicographic order. Can be read using numpy.fromfile(\$Filename\$, dtype=numpy.float32) followed by numpy.reshape()* \
    ├── fluence_maps ...... *contains binary fluence maps of the same number as specified in beamlist. They are all set to zero.* \
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 001.dat\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 002.dat\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── ...\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 020.dat\
    ├── isocenter.txt ...... *text file containing the coordinate of isocenter (in mm)* \
    ├── OAR_target.dat ...... *binary file corresponding to $t$, in the same format as CT.dat. The same below* \
    ├── OAR_weight_correct.dat ...... *with generated new structures such as skin and ring structure, in correspondence with the baseline* \
    ├── OAR_weight.dat ...... *binary file corresponding to $w$* \
    ├── PTV_target.dat ...... *binary file corresponding to $t_0$* \
    ├── PTV_weight.dat ...... *binary file corresponding to $w_0$* \
    ├── shape.txt ...... *text file containing the shape of the CT. For all 6 patients, we use isotropic voxel size of 2mm*
  * patient1_dicom\
    ├── beamlist_annealing.txt ...... *The same as beam_annealing_varianIEC.csv, for baseline input* \
    ├── beamlist_BOO.txt ...... *The same as beam_angles_VarianIEC.csv* \
    ├── beamlist_full.txt ...... *All candidate beams, 1162 beams for each patient. Less for patient1 for memory constraint* \
    ├── config.json ...... *configuration files for baseline dose calculation* \
    ├── dicom ...... *dicom folders* \
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 001.dcm\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 002.dcm\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── ...\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── 197.dcm\
    │ &nbsp;  &nbsp;  &nbsp; &nbsp;├── RTstructure.dcm ...... *contains anatomical structure* \
    ├── structures.json ...... *a list of anatoimcal structures*

### Code execution
  The code is tested on a Ubuntu 18.04.4 LTS, with Intel(R) Core(TM) i9-7900X CPU with 20 cores and 62G memory, and 4 NVIDIA GeForce GTX 1080 ti GPUs, each with 11G memory (although only one is used). CUDA version 10.1.243. Assume the project is under ```~/EndtoEnd```. To compile:
  ```console
  user@host:~/EndtoEnd$ bash configure.sh
  user@host:~/EndtoEnd$ bash build.sh
  ```
  These commands should generate a ```./build``` folder with several executables in it, including ```optimize_annealing_constrain```, ```optimize_stationary_smoothness``` and ```dose_calculation```. If not generated successfully, there might be some package dependency issue, which users are supposed to deal with.

  ```dose_calculation``` takes in beam list and fluence maps of the same number. If the folder of the fluence maps is not specified, it initializes the fluence map with ones. An example bash script is ```./examples/dose_calculation.sh```. Before running the script, the user should revise the arguments. The author suggests using the default fluence map configurations (as is in the file). There are some path arguments such as ```--spectrum-path```, ```ATheta-path``` and so on. They contain the physical parameters of the X-ray. Users can find them in ```./kernels```. Specifically, the user should specify ```--output-folder```, which is where the output files will be stored, and will be created if it doesn't exist. Some other auguments such as ```--iterations``` and ```--step-size``` is meaningless.

  After the parameters are set correctly, to execute the code:
  ```
  user@host:~/EndtoEnd$ bash ./examples/dose_calculation.sh
  ```
  It will generate files ```dose001.dat```, ```dose002.dat``` and so on in the specified output folder. They are dose matrices of size height\*width\*slices'. It should be noted that slice' is the minimum multiple of 8, which is equal or greater than slice. For example, if slice = 132, then slice' = 136. This is to ensure GPU memory alignment. 

  ```optimize_stationary_smoothness``` optimizes the fluence maps with beam angles fixed. The input arguments are almost the same as above. The user should additionally specify ```--iterations``` and ```--step-size``` for optimization. It will generate several files in the specified output folder:
  * ```DoseLoss.dat```: a 32-bit float array of length euqaling to iterations, corresponding to the first two terms in the loss function above.
  * ```fluence_maps```: the folder containing the optimized binary fluence maps. In the same format as the input fluence maps.
  * ```SmoothnessLoss.dat```: a 32-bit float array of length equaling to iterations times number of beams, in lexicographic order.
  * ```totalDose.dat```: a 32-bit float array of the final dose matrix, of size height\*width\*slices'.

  To execute, run
  ```console
  user@host:~/EndtoEnd$ bash ./examples/optimize_stationary_smoothness.sh
  ```

  ```optimize_annealing_constrain``` collectively optimizes the beam orientations and fluence maps, with simulated annealing and gradient descent, respectively. It generates several files in the specified output directory:
  * ```zenith.dat``` and ```azimuth.dat```, each of which is a 32-bit float array of size iterations \* num_beams * 2, for two possibilities per beam per iteration: take or not take.
  * ```PerturbationLoss.dat```: float array of the same size as zenith.dat and azimuth.dat, containing losses of the two possibilites.
  * ```possibility.dat```: probabilities of taking the step. Of shape iterations \* num_beams.
  * ```taken.dat```: bool array indicating whether the step has been taken or not, of shape iterations \* num_beams.
  * ```DoseLoss.dat```, ```fluence_maps```, ```SmoothnessLoss.dat```, ```totalDose.dat```: the same as above.
  
  To execute, run
  ```
  user@host:~/EndtoEnd$ bash ./examples/optimize_annealing_constrain.sh
  ```

[FCBB]: https://pubmed.ncbi.nlm.nih.gov/21081826/