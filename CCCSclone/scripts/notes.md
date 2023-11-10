## This is a note I take for the original CCCS dose calculations
According to my understanding, the `main.cu` file is mainly to prepare for the variables and arguments that are passed to the dose calculation kernel. Here I take a note of the variables it prepares.

The variable `tdata` is a list of `*DEVICE _THREAD_DATA`, which has the following attributes:
* `device_beam_arr`. This attribute contains the list of all beams that are to be processed on the device. Each beam contains information about its angle, SAD, and isocenter.
* `device_nbeams`: the number of beams it contains, or the size of the `device_beam_arr`.
* `nrays`: the number of convolution directions per beamlet. calculated as:
    ```
    int nrays = constants->nphi * (constants->ntheta / 2);
    ```
* `deviceid`: nominally, device id.
* `n_unpackstreams`: Number of threaded workers to use during beamlet sparsify.

Besides the parameters that are passed as arguments depending on the beam index, there are some global arguments.
* `mono_kernels`: it essentially contains information of a list of beams. Each For each beam, there are associated properties: `energy`, `fluence`, `mu`, `mu_en`, and `kernel`. Each kernel maintains 5 categories: `primay`, `first_scatter`, `second_scatter`, `multiple_scatter`, `brem_annih`. And the total value of the beam is the sum of all the categories above.

    It has a method to make a poly-kernel. Specifically, it firstly generates the total kernel per each category, as a weighted sum across all mono-energetic kernels. These kernels are weighted by `fluence * energy * mu`. The total kernel is then a sum from all categories.
* `constants`: it contains other stuff.
    * `size`: dimensions of data in voxels after resampling. (103, 103, 103)
    * `voxel`: size of voxel in mm. (0.25, 0.25, 0.25)
    * `start`: 3D coordinates of data origin
    * `calc_bbox_start`: dose bounding box start indices. (1, 1, 1)
    * `calc_bbox_size`: dose bounding box dimension. (101, 101, 101)
    * `max_rev_size`: REV convolution array dimensions. (800, 800, 800)
    * `rev_latspacing`: Convolution ray lateral spacing. 0.25
    * `rev_longspacing`: convolution step spacing. 0.25
    * `kernel_extent`: dose kernel radius truncate distance (cm). 4.000000
    * `ss_factor`: Terma super-sampling factor. (8)
    * `nphi, ntheta, nradii`: CCK Kernel Dimensions. (8, 8, 24)
    * `penumbra`: margin added to rev in fullbeam calc. (1.000000)
    * `beam_count`: number of beams to pick up from `beam_list.txt`
    * `beam_spectrum`: beam energy spectrum file to use.
    * `target_structure`: Name of selected target contour.

In the method `initCudaConstandTex`, three textures are generated.
* `texKern`. Dimension: `(nradii, ntheta)`, source: `kern.total_matrix`, which is read from file `Paths::Instance()->cum_kernel_file();`.
* `texSpectrum`. Dimension: `(mono->nkernels, 4)`. For eack kernel category, `(fluence, energy, mu_en, mu)`.
* `texDens`. Dimension: `(103, 103, 103)`, source: `datavols->density`.

Dose calculation workflow:
* iterate over beams

    Filter valid beamlets

    Calculate the start coordinates, end coordinates, beamlet angles and beamlet isocenters for each beamlet.

    Use the maximum beamlet size to determine the pillar dimension and size, `hPG.pillarSize`, and `hPG.pillarDims`, respectively, where the `y` dimension is the longitudinal dimension, `x` and `z` dimensions corresponds to the two fluence map dimensions.

    Calculate the pillar start coordinates, `hPG.pillarStartCoords`, for each beamlet.

    Group the beamlets into batches (1 batch, no effect).

    Group the pillars into a near-square configuration, with dimensions `(bPG.numPillars.x, bPG.numPillars.y)`, where `bPG` is the pillar grid for the current batch.

    `bPG.gridDims = (bPG.numPillars.x * bPG.pillarDims.x + bPG.wallThickness,` \
    `bPG.pillarDims.y + bPG.wallThickness,` \
    `bPG.numPillars.y * bPG.pillarDims.z + bPG.wallThickness);`

    `bPG.gridSize = make_float3(bPG.gridDims) * rev_voxelsize`

    `rayBlock = dim3{1, 32, 3};`

    * For each ray

        Retrieve the ray angles `theta` and `phi`.

        <!-- Project the pillar volume (dimension described above) onto the `REV` frame, and calculate the dimensions of the bounding boxes, which are aligned to the coordinate axis of the `REV` frame. `rev[rr].size` is the physical size (cm) of the bounding box, with order `(y, z, x)`, while `(rev[rr].min_coords, rev[rr].max_coords)` are the maximum and minimum coordinates of the bounding box, in physical size. -->

        Project the pillar volume, `bPG.gridSize` to the `REV` frame, giving `rev[rr].min_coords`, `rev[rr].max_coords`, and `rev[rr].size`, where the former two are in order `(x, y, z)`, the latter one in order `(y, z, x)`

        `rayGrid = (rev[rr].size.z, rev[rr].size.x/rayBlock.y, rev[rr].size.y/rayBlock.z)`. Considering the order of `rev`, the dimension of `rayGrid` should be `(x, y, z)`

        `convBlock[rr] = (rev[rr].size.x, 1, 1)`\
        `convGrid[rr] = (1, rev[rr].size.y, rev[rr].size.z)`\
        Considering the order of `rev`, the dimensions these two are both `(y, z, x)` 

    `tileBlock=(32, 4, 1)`, `packedGrid=(bPG.gridDims.y / tileBlock.x, bPG.gridDims.z / tileBlock.y, bPG.gridDims.x / tileBlock.z)`, order: `(y, z, x)`

    `device_data[deviceid].dose` is of size `bPG.pillar_grid_nvoxels()`, which is the product of `bPG.gridDims`.

    * For each ray

        `threadIdx` is arranged in `(z, x, y)` order, in accordance with `rayGrid`.

        Threads are indexed with `(tX, tY, tZ)`, which is of order `(x, y, z)`. So `tY` is the longitudinal direction.

        For each thread, calculate its coordinate in the pillar, `pbev_idx = rotate(rev[rr].min_coords)/rev_voxelsize + (tZ, tX, tY), kernel_theta, kernel_phi) + 0.5`. `(tZ, tX, tY)` is of order `(x, y, z)`

        Calculate the dose in the `REV` frame, using CK algorithm.

        Project the `REV` dose to the `BEV` dose, or pillar dose. \
        Threads order: `(y, z, x)` \
        `(tX, tY, tZ)` order: `(y, z, x)` \
        `(pX, pY, pZ)` order: `(x, y, z)`

        Create a texture memory from the `BEV` dose.

        Create a texture memory from the calculated `BEV` dose. \
        Project the `BEV` dose to the patient volume