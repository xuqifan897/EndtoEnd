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

## This project is still going on and continuously updated.

[FCBB_link]: https://pubmed.ncbi.nlm.nih.gov/21081826/