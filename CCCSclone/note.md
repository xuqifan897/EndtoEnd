This is a note to understand the original CCCS algorithm.

* The variable `hPG` represents beamlet information.
    ```
    hPG.numPillars.x = static_cast<int>(std::ceil(std::sqrt(float(hPG.numBeamlets))));
    hPG.numPillars.y = static_cast<int>(ceil(float(hPG.numBeamlets) / hPG.numPillars.x));
    hPG.gridDims = make_int3(
        hPG.numPillars.x * hPG.pillarDims.x + hPG.wallThickness,
        hPG.pillarDims.y + hPG.wallThickness,
        hPG.numPillars.y * hPG.pillarDims.z + hPG.wallThickness
    );
    hPG.gridSize = make_float3(hPG.gridDims) * rev_voxelsize;
    ```
* `rev[rr]` represents coordinates of the current beam (a set of beamlets) in the view of ray `rr`.
    * `rev->size` is the beam range in the view of ray, unitless, in voxel, order: `yzx`.
    * `rev->min_coords` and `rev->max_coords` are its minimum and maximum coordinates, in physical units, order: `xyz`.
    ```
    rev->size = make_uint3(
        static_cast<unsigned int>( ceil(fabsf(pbev_max.y-pbev_min.y) / constants->rev_longspacing) ),
        static_cast<unsigned int>( ceil(fabsf(pbev_max.z-pbev_min.z) / constants->rev_latspacing) ),
        static_cast<unsigned int>( ceil(fabsf(pbev_max.x-pbev_min.x) / constants->rev_latspacing) )
    );
    rev->min_coords = pbev_min;
    rev->max_coords = pbev_max;
    ```
* `conBlock` is `[rev[rr].size.x + overwrite, 1, 1]`\
    `conGrid` is `[1, rev[rr].size.y + overwrite, rev[rr].size.z + overwrite]`
* `rayGrid[rr]` is to partition `rev[rr].size` to `rayBlock`'s.
    ```
    rayGrid[rr] = make_uint3(rev[rr].size.z,
        static_cast<uint>(std::ceil(static_cast<float>(rev[rr].size.x) / rayBlock.y)),
        static_cast<uint>(std::ceil(static_cast<float>(rev[rr].size.y) / rayBlock.z)));
    ```
* In kernel `cudaBeamletRaytrace`
    * `tX`, `tY`, `tZ` are coordinates inside the ray view. They are in the order `yzx`.
    * `pbev_idx` is the coordinates in the beam view, in the order `xyz`.
    * `pillar_idx` is the index of the pillar (beamlet), `pidx` is the beamlet index.
        ```
        int2 pillar_idx = make_int2(
            __float2int_rd(__fdiv_rn(pbev_idx.x, f_pg_pillarDims.x)),
            __float2int_rd(__fdiv_rn(pbev_idx.z, f_pg_pillarDims.z))
        );
        pidx = pillar_idx.x + pillar_idx.y * pg_numPillars.x;
        ```
    * `in_pillar_idx` is the coordinate inside the pillar, unitless.
        ```
        in_pillar_idx = make_float3(
            fmodf(pbev_idx.x, f_pg_pillarDims.x),
            fmodf(pbev_idx.y, f_pg_pillarDims.y),
            fmodf(pbev_idx.z, f_pg_pillarDims.z)
        );
        ```
    * `g_coords` is the coordinates of the current ray voxel in the patient frame of view, in physical length.
        ```
        g_coords = d_inverseRotateBeamAtOriginRHS(
            d_inverseRotateBeamletAtOriginRHS(sub_in_pillar_idx, beamletAngles.x, beamletAngles.y),
            beam_azimuth, beam_zenith, beam_coll);
        g_coords = __fmaf_rn(g_coords, rev_voxelsize, pillarStartCoords);
        ```
    * This function initializes the density and Terma for the current ray in the ray view. `out_idx` is in the order `(tZ, tY, tX)`, or `xzy`
* In kernel `PackRowConvolve`
    * `(pack_tX, pack_tY, pack_tZ) = [threadIdx.x, blockIdx.y, blockIdx.z]`, in the order `(rev.size.x, rev.size.y, rev.size.z)`, or `yzx`.
    * `mem_idx` in the order `(pack_tZ, pack_tY, pack_tX)`, or `xzy`.
    * This function calculates the dose in the ray view, and store the dose at `surfDoseObj`, at coordinate `(pack_tX, pack_tY, pack_tZ)`, or `yzx` (from least significant to most significant).

* `packedGrid` is to partition the whole beam into `tileBlock`'s. `tileBlock = [32, 4, 1]`.
    ```
    dim3 packedGrid = dim3(
        static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.y) / tileBlock.x)),
        static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.z) / tileBlock.y)),
        static_cast<uint>(std::ceil(static_cast<float>(hPG.gridDims.x) / tileBlock.z))
    );
    ```
* In kernel `PackedREVtoBEVdose`
    * `(tX, tY, tZ)` is in the order `(hPG.gridDims.y, hPG.gridDims.z, hPG.gridDims.x)`, or `yzx`.
    * `(pX, pY, pZ)` is in the order `(tZ, tX, tY)`, or `xyz`.
    * `rev_idx` is to transform the current coordinates from the beam view to the ray view.
        ```
        float3 rev_idx = __fmaf_rn(-revStart, __frcp_rn(rev_voxelsize), 
            d_rotateKernelAtOriginRHS(make_float3(pX, pY, pZ), kern_theta, kern_phi));
        ```
        `rev_idx` is in the order `(pX, pY, pZ)`, or `xyz`
    * `dose_coeff` is read from the texture `texPackedREVDose`, with coordinates `(rev_idx.y + 0.5f, rev_idx.z + 0.5f, rev_idx.x + 0.5f)`. The coordinates are in the order `yzx`, which is consistent with the order of `(pack_tX, pack_tY, pack_tZ)` in the kernel `PackRowConvolve`.
    * `dose_coeff` is added to the beam dose, at index `(tZ, tY, tX)`, or `xzy` (from most significant to least significant).