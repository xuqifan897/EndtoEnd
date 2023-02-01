# Dose Calculation (GPU)
This project aims to be an efficient gpu-based External-Beam Radiation Therapy (EBRT) dose calculation software capable of producing both _beam_- and _beamlet-based_ dose volume data. 

## Motivation
While _beam-based_ deterministic dose calculation is mostly considered a "solved" problem (since many clinical software packages can efficiently do this in little time with little compute resources), the problem of large-scale _beamlet-based_ dose calculation is still critically unhandled. The primary contribution of this project is to improve upon existing GPU _beamlet-based_ deterministic dose calculation methods using a novel _context-based_ dose approximation, while also enabling additional parallelization across a scalable network of multi-GPU compute nodes. The result is considerable improvements to _beamlet-based_ dose calculation speeds, which enable consideration of hundreds or even thousands of External X-ray beams when optimizing treatment plans. In most cases, runtimes of under 25 minutes for ~1000 beams can be expected with just a single node and 4x GPUs. Even greater speed gains are achievable with increased GPU hardware availability.

## History
This package is based on a _beam-based-only_ GPU implementation of Jack Neylon's (UCLA Radiation Oncology) Non-voxel Based (NVB) extension ([paper][neylon-nvb-dosecalc-paper]) of the Pinnacle collapsed cone convolution superposition (CCCS) code for computing x-ray dose in external beam radiotherapy. The project is currently maintained by Ryan Neph (PhD Student - UCLA Radiation Oncology). Initial implementation of beamlet-based modifications to the NVB code was carried out by Ryan Neph with help from Cheng Ouyang (Imperial College London) during his summer internship at UCLA in 2017. Modifications have been made to produce sparsified beamlet-specific dose distributions for use in an in-house inverse treatment planning optimization software (Dr. Sheng's 4pi method, among others). Details on the conversion to a beamlet-based algorithm are provided in the recent [paper][neph-gpu-dosecalc-paper].

## Getting Started
Before running the dose calculation, one must first run a preprocessing step which configures the settings and prepares the structured inputs for the eventual dose calculation:

* __`dosecalc-preprocess`__ - handles import of Dicom data and manipulation of algorithm inputs into standard binary formats for use by `dosecalc-beamlet`. This will be crucial for setting up task data for separate networked compute nodes.

After preprocessing is complete, one can select between a _beamlet-based_ or _beam-based_ dose calculation and run the corresponding program. All inputs to the dose calculation program are assumed to be located in the right place on the filesystem after running the preprocessing step, so the only arguments that can be provided to the following programs are specific to each one and generally allow changes to the outputs rather than the calclulation procedure itself (which should be specified during preprocessing instead):

* __`dosecalc-beamlet`__ - reads binary task data from `dosecalc-preprocess` and runs multi-GPU beamlet-based dose calculation for constructing the "Dose Coefficient Matrix". Data is packaged into a single [hdf5 binary result file](https://gitlab.com/shenglab/dose-calculation/wikis/Beamlet-DoseCalc-Data-Output-Specification). 

    _\- or -_

* __`dosecalc-beam`__ - reads binary problem data from `dosecalc-preprocess` and schedules gpu kernel execution for efficient dose calculation and accumulation. This is different from `dosecalc-beamlet` in that it produces a single dose volume from all beams.

## Installation
By far the easiest way to begin using this program is via Docker. For convenience, a Dockerfile is provided which defines a fully-working build of a local Docker Image. A pre-built image is also provided on the gitlab registry for the project (<https://gitlab.com/shenglab/dose-calculation/container_registry>) and can be pulled using the command: 

    docker pull registry.gitlab.com/shenglab/dose-calculation:stable

### Dependencies
#### Software
* [CMAKE] >=v3.14
* [CUDA Toolkit] >=v7.5 (support of c++11 standard), and supported nvidia graphics drivers (find compatibility [here](https://docs.nvidia.com/deploy/cuda-compatibility/index.html#binary-compatibility))
* DCMTK [>=v3.6.2][dcmtk] / [latest dev-snapshot][dcmtk-dev] (be mindful of GCC version compatability; latest dev snapshot should be okay with GCC \>=4.8.x)
* [Boost C++ Library] >=1.30 (Ubuntu: `sudo apt-get install libboost-all-dev`)

** Note, by using the Docker image, the software dependencies above are already taken care of but the following additional dependencies are required. (These may already be configured on your system if you have used Docker for GPU-based execution before):
* [Docker] >= 1.12 (latest is recommended)
* [nvidia-docker] v2.0

#### Hardware: 
* 1\+ CUDA-capable NVIDIA GPU(s) [(compute capability >=3.0 / >=Kepler)][cuda-cc]
* 2GB\+ GPU Memory Recommended

### Tested Platforms
* Ubuntu 14.04 LTS (kernel 3.19.x) / GCC 4.9.x / DCMTK development snapshot 3.6.1-20170228
* Ubuntu 16.04 LTS (kernel 4.4.x) / GCC 5.4.x / DCMTK development snapshot 3.6.1-20170228
* Ubuntu 16.04 LTS (kernel 4.4.x) / GCC 5.4.x / DCMTK 3.6.2

[neph-gpu-dosecalc-paper]: https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.13651
[neylon-nvb-dosecalc-paper]: https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4895822

[CUDA Toolkit]: https://developer.nvidia.com/cuda-toolkit-archive
[dcmtk]: http://dicom.offis.de/dcmtk.php.en
[dcmtk-dev]: http://dicom.offis.de/download/dcmtk/snapshot/
[CMAKE]: https://cmake.org/download/
[cuda-cc]: https://developer.nvidia.com/cuda-gpus
[Docker]: https://docs.docker.com/install/
[nvidia-docker]: https://github.com/nvidia/nvidia-docker/wiki/Installation-(version-2.0)
[Boost C++ Library]: https://www.boost.org/
