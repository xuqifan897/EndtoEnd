#include "nvbbRayConvolve_device.cuh"
#include "DoseCalcAlgo/cudaSiddon.cuh"

#define INVALID     -1e-12;
#define subvoxidx   (threadIdx.x)
#define subvoxcount (blockDim.x)
#define s_index(ii) threadIdx.y+blockDim.y*(threadIdx.z+blockDim.z*ii)
#define SUBTERMA s_subterma[s_index(subvoxidx)]

#define SET_INACTIVE_THREAD { \
    s_subterma[s_index(0)]=-1; \
    valDens = INVALID; \
    valTerm = 0.f; \
}
#define CHECK_ACTIVE_THREAD (s_subterma[s_index(0)]>=0)
// kernel to trace through density and terma volumes
// resampling along convolution direction into packed pillar grid configuration
#define raytrace_write(arr, val) arr[out_idx] = val;
__global__ void
cudaBeamletRaytrace(
        float               *packbDens,                 // storage for BEV resampled density
        float               *packbTerm,                 // storage for BEV resampled terma
        float3              beam_source,
        float2              beamlet_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        float3              f_pg_gridDims,
        float3              f_pg_pillarDims,
        int                 pg_numBeamlets,
        int                 pg_wallThickness,
        int2                pg_numPillars,
        int*                pg_beamletIdx,
        float2*             pg_beamletAngles,
        float3*             pg_pillarStartCoords,
        float3*             pg_beamletIsocenters,
        float               kern_theta,
        float               kern_phi,                   // kernel rotation angles
        float3              revStart,                   // REV volume limit coords in XYZ coord system
        uint3               rev_size,                   // size of resampled REV data volumes
        uint3               max_rev_size,
        float3              start,
        float3              voxel,
        float3              rev_voxelsize,
        float3              f_calc_bbox_start,
        cudaTextureObject_t tex3Ddens,
        cudaTextureObject_t tex3Dter,

        float* d_fluence_map,
        float3 f_size,
        uint3 calc_bbox_size,
        float beamhard_correct,
        float3 beam_direction,
        float3 beam_isocenter,
        float beam_sad,
        uint2 fmap_size,
        uint dkerns,
        cudaTextureObject_t specTex,
        int requested_beamletnum = -1                   // if <0, raytrace for all beamlets
) {
    // supersampling setup
    extern __shared__ float s_subterma[];
    SUBTERMA = 0.f;
    __syncthreads();

    //------------------------------------

    float valDens = 0.f;
    float valTerm = 0.f;

    int out_idx;
    int pidx, beamletnum;
    float3 in_pillar_idx;
    if (subvoxidx == 0) {
        // thread launched for each voxel in REV space
        // Pay attention: we write the data along the x-axis but this maps to the y-axis in the pBEV space
        // some uses of the axes will thus be transposed
        {
            int tX = threadIdx.y + blockIdx.y * blockDim.y;
            int tY = threadIdx.z + blockIdx.z * blockDim.z;
            int tZ = blockIdx.x;
            /* // pBEV indices */
            /* int pX = tZ; */
            /* int pY = tX; */
            /* int pZ = tY; */
            // check for out of REV bounds
            if (tX >= rev_size.x || tY >= rev_size.y || tZ >= rev_size.z ) { SET_INACTIVE_THREAD; }
            out_idx = tX + (max_rev_size.x * (tY + (max_rev_size.y * tZ)));

            // convert REV indices to pBEV indices, shift to 0 at edge of voxel for position tests
            float3 pbev_idx = inverseRotateKernelAtOriginRHS(__fmaf_rn(revStart, __frcp_rn(rev_voxelsize), make_float3(tZ, tX, tY)), kern_theta, kern_phi)+0.5f;

            // check for out of pBEV bounds
            if (pbev_idx.x < 0.f || pbev_idx.y < 0.f || pbev_idx.z < 0.f
                || pbev_idx.x >= f_pg_gridDims.x
                || pbev_idx.y >= f_pg_gridDims.y
                || pbev_idx.z >= f_pg_gridDims.z )
            {
                /* raytrace_write(packbDens, -3.f); */
                SET_INACTIVE_THREAD;
            }
            // get pillar index
            if (CHECK_ACTIVE_THREAD){
                int2 pillar_idx = make_int2(
                        __float2int_rd(__fdiv_rn(pbev_idx.x, f_pg_pillarDims.x)),
                        __float2int_rd(__fdiv_rn(pbev_idx.z, f_pg_pillarDims.z))
                        );
                pidx = pillar_idx.x + pillar_idx.y * pg_numPillars.x;

                // check if inactive pillar
                if (pidx >= pg_numBeamlets) {
                    /* raytrace_write(packbDens, -2.f); */
                    SET_INACTIVE_THREAD;
                }
            }

            // only extract context for current beamlet
            if (CHECK_ACTIVE_THREAD) {
                beamletnum = pg_beamletIdx[pidx];
                if (requested_beamletnum >=0 && beamletnum != requested_beamletnum) {
                    SET_INACTIVE_THREAD;
                }
            }

            if (CHECK_ACTIVE_THREAD) {
                // check for wall
                in_pillar_idx = make_float3(
                        fmodf(pbev_idx.x, f_pg_pillarDims.x),
                        fmodf(pbev_idx.y, f_pg_pillarDims.y),
                        fmodf(pbev_idx.z, f_pg_pillarDims.z)
                        );
                if (in_pillar_idx.x < pg_wallThickness || in_pillar_idx.y < pg_wallThickness || in_pillar_idx.z < pg_wallThickness) {
                    /* raytrace_write(packbDens, -1.f); */
                    SET_INACTIVE_THREAD;
                }
            }
        }
    }
    __syncthreads();

    // Only proceed if thread is active
    if (CHECK_ACTIVE_THREAD) {
        // Fetch density and terma from textures for remaining threads;
        // pillarStartCoords provide anchor to RCS and center for beamlet divergence and beam angle rotations
        // TODO: move pillarStartCoords loading into shared memory
        // TODO: move beamletAngles loading into shared memory

        in_pillar_idx -= 0.5f; // shift back to 0 at center of voxel

        float3 g_coords;
        float2 beamletAngles = pg_beamletAngles[pidx];
        float3 pillarStartCoords = pg_pillarStartCoords[pidx];
        float3 sub_in_pillar_idx = in_pillar_idx;

        // compute anti-aliased terma
        // get sub-vox center position
        if (subvoxcount > 1) {
            int ss_factor = cbrtf(subvoxcount);
            int sub_X = subvoxidx % ss_factor;
            int sub_Y = (subvoxidx / ss_factor)%ss_factor;
            int sub_Z = subvoxidx / (ss_factor*ss_factor);

            float3 subvoxsize = make_float3(1.f/(float)ss_factor);
            sub_in_pillar_idx = make_float3(
                    in_pillar_idx.x - 0.5f + (sub_X+0.5f)*subvoxsize.x,
                    in_pillar_idx.y - 0.5f + (sub_Y+0.5f)*subvoxsize.y,
                    in_pillar_idx.z - 0.5f + (sub_Z+0.5f)*subvoxsize.z );
        }

        g_coords = inverseRotateBeamAtOriginRHS(
                inverseRotateBeamletAtOriginRHS(sub_in_pillar_idx, beamletAngles.x, beamletAngles.y),
                beam_azimuth, beam_zenith, beam_coll);
        g_coords = __fmaf_rn(g_coords, rev_voxelsize, pillarStartCoords);

        // XXX THIS VERSION CAUSES incorrect context alignment (wrong beamlet divergence correction)
        /* g_coords = inverseRotateBeamletAtOriginRHS( */
        /*                    inverseRotateBeamAtOriginRHS(sub_in_pillar_idx*rev_voxelsize, beam_azimuth, beam_zenith, beam_coll), */
        /*                beamletAngles.x, beamletAngles.y) + pillarStartCoords; */

        /* float3 g_term_indices = float3{ */
        /*     __fsub_rn(__fdiv_rn(__fsub_rn(g_coords.x, start.x), voxel.x), f_calc_bbox_start.x), */
        /*     __fsub_rn(__fdiv_rn(__fsub_rn(g_coords.y, start.y), voxel.y), f_calc_bbox_start.y), */
        /*     __fsub_rn(__fdiv_rn(__fsub_rn(g_coords.z, start.z), voxel.z), f_calc_bbox_start.z) }; */
        /* valTerm = tex3D<float>( tex3Dter,  g_term_indices.x, g_term_indices.y, g_term_indices.z); */

        // Split into subvoxels for supersampled terma averaging
        float3 end;
        end.x = start.x + f_size.x * voxel.x;
        end.y = start.y + f_size.y * voxel.y;
        end.z = start.z + f_size.z * voxel.z;
        SUBTERMA = cudaSiddon_device(
                g_coords,
                d_fluence_map,
                start,
                end,
                voxel,
                beamhard_correct,
                beam_source,
                beam_direction,
                beam_isocenter,
                beam_sad,
                beamlet_size,
                fmap_size,
                beam_azimuth,
                beam_zenith,
                beam_coll,
                dkerns,
                tex3Ddens,
                specTex,
                beamletnum
                );

        if (subvoxidx == 0) {
            // compute density at center of voxel
            g_coords = inverseRotateBeamAtOriginRHS(
                    inverseRotateBeamletAtOriginRHS(in_pillar_idx*rev_voxelsize, beamletAngles.x, beamletAngles.y),
                    beam_azimuth, beam_zenith, beam_coll) + pillarStartCoords;
            /* g_coords = inverseRotateBeamAtOriginRHS( */
            /*         inverseRotateBeamletAtOriginRHS(in_pillar_idx, beamletAngles.x, beamletAngles.y), */
            /*         beam_azimuth, beam_zenith, beam_coll); */
            /* g_coords = __fmaf_rn(g_coords, rev_voxelsize, pillarStartCoords); */
            float3 g_dens_indices = __fdiv_rn(__fsub_rn(g_coords, start), voxel);
            valDens = tex3D<float>( tex3Ddens, g_dens_indices.x+0.5f, g_dens_indices.y+0.5f, g_dens_indices.z+0.5f );
        }
    }

    // ---------------------------------
    __syncthreads();

    if (subvoxidx==0) {
        if (CHECK_ACTIVE_THREAD) {
            float avgterma = 0.f;
            for (int ii=0; ii<subvoxcount; ii++) {
                avgterma += s_subterma[s_index(ii)];
            }
            valTerm = avgterma / subvoxcount;
        }

        raytrace_write(packbDens, valDens);
        raytrace_write(packbTerm, valTerm);
    }
}


// perform a line convolution in both directions along each ray
// each block is parallel ray and each thread is sample point along ray (spaced by delR)
// beamlet version
#define pack_tX threadIdx.x
#define pack_tY blockIdx.y
#define pack_tZ blockIdx.z
#define SHM_DENS(v) shmBlock[v]
#define SHM_TERM(v) shmBlock[rev_size.x + v]
__global__ void
PackRowConvolve( float *bDens,
             float *bTerma,
             cudaSurfaceObject_t surfDoseObj,
             float f_kern_wt_idx,
             uint3 rev_size,                     // size of resampled REV data volumes
             uint3 max_rev_size,
             float delr,
             uint nradii,
             uint ntheta,
             uint nphi,
             cudaTextureObject_t kernTex,
             bool TERMA_ONLY
) {
	// initialize doseArray to 0 since last convolution
    surf3Dwrite(0.f, surfDoseObj, sizeof(float) * pack_tX, pack_tY, pack_tZ, cudaBoundaryModeZero);

	// each row of data now corresponds to a single ray, as was sampled in the cudaRaytrace kernel
	// a block of threads is launched for each parallel REV ray with a thread for each sampling point along the ray
    int mem_idx = pack_tX + (max_rev_size.x * (pack_tY + (max_rev_size.y * pack_tZ)));

	// dynamically allocated shared memory
    extern __shared__ float shmBlock[];

    // Move density and terma data into shared memory space on each block
    if (pack_tX < rev_size.x && pack_tY < rev_size.y && pack_tZ < rev_size.z ) {
        SHM_DENS(pack_tX) = bDens[mem_idx];
        SHM_TERM(pack_tX) = bTerma[mem_idx];
    }
    __syncthreads();

    // check for out of REV bounds
    if (pack_tX >= rev_size.x || pack_tY >= rev_size.y || pack_tZ >= rev_size.z ) { return; }

    // check for pillar wall, invalid pillar, outside of pillar
    if (SHM_DENS(pack_tX) < 0.0f) { return; }

    if (TERMA_ONLY) {
        surf3Dwrite(SHM_TERM(pack_tX), surfDoseObj, sizeof(float) * pack_tX, pack_tY, pack_tZ, cudaBoundaryModeZero);
        return;
    }

	// determine convolution direction index, for sampling the poly-energetic dose deposition kernel
    float P = f_kern_wt_idx;
    float recip_P = ntheta - f_kern_wt_idx - 1;

	// initialize variables
    float sum = 0.f;
    int kk = pack_tX;
    /* float r_eff = __fmul_rn(0.5f, __fmul_rn(delr, SHM_DENS(kk))); // start at first voxel wall */
    float r_eff = 0.f;
    int idx_rad = 0; // (far bound)
    float last_cumval = 0.f;

	// accumulate convolution contributions in one direction along the row of data
    // This the "forward" direction, such that an interaction point deposits dose at this thread's assignment
    // voxel by traveling in the forward kernel direction
    do {
		// break loop early if sampling outside active calculation volume
        if (kk < 0) { break; }
        /* if ( r_eff >  devConstants.kernel_extent ) { break; } */
        if ( SHM_DENS(kk) < 0.f ) { break; } // stop convolution along this direction when boundary is reached

        // Find effective radiological distance to voxel kk (heterogeneity correction)
        float delr_eff = delr*SHM_DENS(kk);
        r_eff += 0.5f*delr_eff; // move to center of this voxel (interaction at center)

        if (r_eff > KERN_RADII[nradii - 1]) { break; }

        // find bounding kernel radii for interp. (high wall)
        while (KERN_RADII[idx_rad] < r_eff) {
            idx_rad++;
        }
        /* float vr = 0.f; */
        /* for (int ll=0; ll<nradii; ll++) { */
        /*     if (ll == 0) { */
        /*         vr = 0.5*KERN_RADII[0]; */
        /*     } else { */
        /*         vr = 0.5*(KERN_RADII[ll] + KERN_RADII[ll-1]); */
        /*     } */
        /*     if (vr >= r_eff) { */
        /*         idx_rad = ll; */
        /*         break; */
        /*     } */
        /* } */

        float cumval;
        if (idx_rad <= 0) {
            float interp = __fmul_rn(r_eff, __frcp_rn(KERN_RADII[0]));
            cumval = __fmul_rn(tex2D<float>( kernTex, 0.5f, P + 0.5f), interp);
        } else {
            float interp = __fmul_rn(__fsub_rn(r_eff, KERN_RADII[idx_rad-1]), __frcp_rn(__fsub_rn(KERN_RADII[idx_rad], KERN_RADII[idx_rad-1])));
            cumval = tex2D<float>( kernTex, idx_rad - 0.5f + interp, P + 0.5f);
        }
        /* cumval = tex2D<float>(kernTex, __int2float_rn(idx_rad) + 0.5f, P + 0.5f); */

        // CCK dose calculation
        sum += (cumval-last_cumval)*SHM_TERM(kk);
        last_cumval = cumval;
        if (kk != pack_tX) {
            r_eff += 0.5f*delr_eff; // move to end of this voxel for next loop
        }
        --kk;
    } while (idx_rad < nradii - 1);

	// reset variables
    kk = pack_tX;
    /* r_eff = __fmul_rn(0.5f, __fmul_rn(delr, SHM_DENS(kk))); */
    r_eff = 0.f;
    idx_rad = 0;
    last_cumval = 0.f;

    // perform the same operation in the opposite direction along the row of data
    // This the "reverse" direction, such that an interaction point deposits dose at this thread's assignment
    // voxel by traveling in the backward kernel direction
    do {
        if (kk >= rev_size.x) { break; }
        /* if ( r_eff > devConstants.kernel_extent ) { break; } */
        if ( SHM_DENS(kk) < 0.f ) { break; }

        float delr_eff = __fmul_rn(delr, SHM_DENS(kk));
        r_eff += 0.5f*delr_eff;

        if (r_eff > KERN_RADII[nradii - 1]) { break; }

        // find bounding kernel radii for interp. (high wall)
        while (KERN_RADII[idx_rad] < r_eff) {
            idx_rad++;
        }
        /* float vr = 0.f; */
        /* for (int ll=0; ll<nradii; ll++) { */
        /*     if (ll == 0) { */
        /*         vr = 0.5*KERN_RADII[0]; */
        /*     } else { */
        /*         vr = 0.5*(KERN_RADII[ll] + KERN_RADII[ll-1]); */
        /*     } */
        /*     if (vr >= r_eff) { */
        /*         idx_rad = ll; */
        /*         break; */
        /*     } */
        /* } */

        float cumval;
        if (idx_rad <= 0) {
            float interp = __fmul_rn(r_eff, __frcp_rn(KERN_RADII[0]));
            cumval = __fmul_rn(tex2D<float>( kernTex, 0.5f, recip_P + 0.5f), interp);
        } else {
            float interp = __fmul_rn(__fsub_rn(r_eff, KERN_RADII[idx_rad-1]), __frcp_rn(__fsub_rn(KERN_RADII[idx_rad], KERN_RADII[idx_rad-1])));
            cumval = tex2D<float>( kernTex, idx_rad - 0.5f + interp,  recip_P + 0.5f);
        }
        /* cumval = tex2D<float>(kernTex, __int2float_rn(idx_rad) + 0.5f, recip_P + 0.5f); */

        sum += (cumval-last_cumval)*SHM_TERM(kk);
        last_cumval = cumval;
        if (kk != pack_tX) {
            r_eff += 0.5f*delr_eff; // move to end of this voxel for next loop
        }
        ++kk;
    } while (idx_rad < nradii - 1);

    // normalize by nphi
    sum /= nphi;

	// write the summed result out using a  object
    surf3Dwrite(sum, surfDoseObj, sizeof(float) * pack_tX, pack_tY, pack_tZ, cudaBoundaryModeTrap);
}


// accumlate convolution results in pBEV grid (appending sum is important here)
#define packedbev2rev_append(val) bev_packed_dose[tX + (pg_gridDims.y * (tY + (pg_gridDims.z * tZ)))] += val;
// pbev indices
#define pX  tZ
#define pY  tX
#define pZ  tY
__global__ void
PackedREVtoBEVdose(
        float*              bev_packed_dose,            // beamlet-packed dose array in BEV orientation
        cudaTextureObject_t texPackedREVDose,           // packed dose array embedded in REV bounding box
        float               kern_theta, float kern_phi, // convolution direction
        float3              revStart,                    // REV volume limit coords in XYZ coord system
        float3              rev_voxelsize,
        int3                pg_gridDims
) {
    // thread launched for every voxel in packed BEV array
    int tX = threadIdx.x + blockIdx.x * blockDim.x;
    int tY = threadIdx.y + blockIdx.y * blockDim.y;
    int tZ = blockIdx.z;

    // check for out of pBEV bounds
    if (pX >= pg_gridDims.x || pY >= pg_gridDims.y || pZ >= pg_gridDims.z) { return; }

    // get packed rev indices (float)
    float3 rev_idx = __fmaf_rn(-revStart, __frcp_rn(rev_voxelsize), rotateKernelAtOriginRHS(make_float3(pX, pY, pZ), kern_theta, kern_phi));
    float dose_coeff = tex3D<float>(texPackedREVDose, rev_idx.y + 0.5f, rev_idx.z + 0.5f, rev_idx.x + 0.5f);
    if (dose_coeff <= 0.f) { return; }

    // sum-write value back to BEV xyz in coalesced manner
    // coalesced global memory write
    packedbev2rev_append(dose_coeff);
}


// unpack a single pillar from the packed BEV dose coefficient array
// output will be XYZ volume bounded by calc_bbox that contains dose contributions from this beamlet (pillar)
// in the standard dicom shape.
#define unpack_write(val) unpacked_dose[tX + (calc_bbox_size.x * (tY + (calc_bbox_size.y * tZ)))] = val;
__global__ void
UnpackBEVDosePillar (
        float*              unpacked_dose,    // output array with size matching calc_bbox for unpacked coeff. from one beamlet
        cudaTextureObject_t texPackedBEVDose, // packed dose array oriented in BEV
        float               beam_sad,
        float3              beam_source,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        float3              start,
        float3              voxel,
        float3              f_calc_bbox_start,
        uint3               calc_bbox_size,
        float3              rev_voxelsize,
        int                 pillarIndex,     // indicates which pillar to unpack
        int                 pillarIndex_x,     // indicates which pillar to unpack
        int                 pillarIndex_y,     // indicates which pillar to unpack
        int3                pg_pillarDims,
        int                 pg_wallThickness,
        int                 pg_pillarBuffer,
        float3              pillarStartCoords, // RCS coords indicating the position of the first pillar voxel
        float2              beamletAngles // beamlet divergence angles attached to this pillar
) {
    // thread launched for every voxel in calc_bbox
    // Get RCS indices (row major sequencing)
    int tX = threadIdx.x + blockIdx.x * blockDim.x;
    int tY = threadIdx.y + blockIdx.y * blockDim.y;
    int tZ = blockIdx.z;

    // check for out of bbox bounds
    if (tX < 0 || tY < 0 || tZ < 0 || tX >= calc_bbox_size.x || tY >= calc_bbox_size.y || tZ >= calc_bbox_size.z ) { return; }

    // bbox(RCS) indices to pbev(intra-pillar) indices
    float3 g_coords = (make_float3(tX, tY, tZ) + f_calc_bbox_start) * voxel + start;
    /* float3 rot_anchor = beam_source-pillarStartCoords; */
    /* float3 p_idx = rotateBeamAtOriginRHS( */
    /*                    rotateBeamletAtOriginRHS(g_coords-pillarStartCoords, beamletAngles.x, beamletAngles.y), */
    /*                    beam_azimuth, beam_zenith, beam_coll */
    /*                )/rev_voxelsize + 0.5f; */
    float3 p_idx = rotateBeamletAtOriginRHS(
                       rotateBeamAtOriginRHS(g_coords-pillarStartCoords, beam_azimuth, beam_zenith,beam_coll),
                       beamletAngles.x, beamletAngles.y
                   )/rev_voxelsize + 0.5f;

    // check for out of pillar or in wall
    if (p_idx.x < (pg_wallThickness+pg_pillarBuffer) || p_idx.y < (pg_wallThickness+pg_pillarBuffer) || p_idx.z < (pg_wallThickness+pg_pillarBuffer)
            || p_idx.x>=(pg_pillarDims.x-pg_pillarBuffer) || p_idx.y>=(pg_pillarDims.y-pg_pillarBuffer) || p_idx.z>=(pg_pillarDims.z-pg_pillarBuffer)) { return; }
    // convert indices from pillar start to indices from array start
    p_idx = make_float3(
            p_idx.x + pillarIndex_x*pg_pillarDims.x,
            p_idx.y,
            p_idx.z + pillarIndex_y*pg_pillarDims.z
            );

    // fetch dose coefficient (remember, pbev axes are transposed)
    float dose_coeff = tex3D<float>(texPackedBEVDose,
            p_idx.y,
            p_idx.z,
            p_idx.x );
    if (dose_coeff <= 0.f) { return; }

    unpack_write(dose_coeff);
}


