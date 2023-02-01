#include "cudaSiddon.cuh"
#include <cassert>

// kernel combines the radiological effective depth calculation (based on Siddon's algorithm)
// with the TERMA calculation
//

// define to make terma add to existing terma array, rather than overwrite all voxels
#define siddon_lidx ( ray_X + calc_bbox_size.x*(ray_Y + calc_bbox_size.y*ray_Z) )
#if 0
#define siddon_write(val) terma[siddon_lidx]+=val
#else
#define siddon_write(val) terma[siddon_lidx]=val
#endif

#define subvoxidx   (threadIdx.x)
#define subvoxcount (blockDim.x)
#define s_index(ii) threadIdx.y+blockDim.y*(threadIdx.z+blockDim.z*ii)

__device__ float cudaRaytrace_device(
        float3              this_voxel,
        float3              start,
        float3              end,
        float3              voxel,
        float3              beam_source,
        cudaTextureObject_t tex3Dvol
) {
    // Siddon's algorithm normalizes the distance between points 1 and 2
    // alpha is the current position along that normalized vector as a scalar amount
    // alpha_min and alpha_max set the bounds of intersection between this vector and the volume of interest
    // limits of alpha parameter - aligned to plane positions not voxel center

    float3 diff = this_voxel-beam_source;
    float3 startdiff = start - beam_source;
    float3 enddiff = end - beam_source;

    // get first and last planar parameterizations
    float a_min, a_max;
    {
        float3 a_first, a_last;
        a_first.x = ((startdiff.x) - 0.5f*voxel.x) / (diff.x);
        a_first.y = ((startdiff.y) - 0.5f*voxel.y) / (diff.y);
        a_first.z = ((startdiff.z) - 0.5f*voxel.z) / (diff.z);
        a_last.x  = ((enddiff.x)   + 0.5f*voxel.x) / (diff.x);
        a_last.y  = ((enddiff.y)   + 0.5f*voxel.y) / (diff.y);
        a_last.z  = ((enddiff.z)   + 0.5f*voxel.z) / (diff.z);

        float3 alpha_min, alpha_max;
        alpha_min.x = fminf(a_first.x, a_last.x);
        alpha_max.x = fmaxf(a_first.x, a_last.x);
        alpha_min.y = fminf(a_first.y, a_last.y);
        alpha_max.y = fmaxf(a_first.y, a_last.y);
        alpha_min.z = fminf(a_first.z, a_last.z);
        alpha_max.z = fmaxf(a_first.z, a_last.z);

        // a_min is alpha of first plane intersection, a_max is location of this threads voxel
        a_min = fmaxf(0,fmaxf(fmaxf(alpha_min.x, alpha_min.y), alpha_min.z));
        a_max = fminf(1,fminf(fminf(alpha_max.x, alpha_max.y), alpha_max.z));
    }

    float rpl = 0.0f;			//radiological path length in cm
    if (!(a_min >= a_max)) {
        float d12 = length(diff);   //distance between ray end points
        float step_x = fabsf(voxel.x/diff.x);
        float step_y = fabsf(voxel.y/diff.y);
        float step_z = fabsf(voxel.z/diff.z);

        // step along the vector, sampling each time we cross any plane (meaning we've entered a new voxel)
        float alpha = a_min;
        float alpha_x = a_min;
        float alpha_y = a_min;
        float alpha_z = a_min;
        float nextalpha;
        float alpha_mid;
        int max_iters = 5000;
        int iter = 0;
        while (alpha < a_max && ++iter < max_iters) {
            // find next intersection plane
            bool valid_x, valid_y, valid_z;
            if (alpha_x >= 0.0f && alpha_x < 1.0f) { valid_x = true; } else { valid_x = false; };
            if (alpha_y >= 0.0f && alpha_y < 1.0f) { valid_y = true; } else { valid_y = false; };
            if (alpha_z >= 0.0f && alpha_z < 1.0f) { valid_z = true; } else { valid_z = false; };
            if (!(valid_x || valid_y || valid_z)) { break; }
            if (valid_x && (!valid_y || alpha_x <= alpha_y) && (!valid_z || alpha_x <= alpha_z)) {
                nextalpha = alpha_x;
                alpha_x += step_x;
            }
            if (valid_y && (!valid_x || alpha_y <= alpha_x) && (!valid_z || alpha_y <= alpha_z)) {
                nextalpha = alpha_y;
                alpha_y += step_y;
            }
            if (valid_z && (!valid_x || alpha_z <= alpha_x) && (!valid_y || alpha_z <= alpha_y)) {
                nextalpha = alpha_z;
                alpha_z += step_z;
            }

            // the total intersection length of previous voxel
            float intersection = fabsf(d12*(nextalpha-alpha));	//intersection is voxel intersection length

            if (intersection>=0.001f) { //do not process unless > 0.01 cm
                alpha_mid = (nextalpha + alpha)*0.5f;  //midpoint between intersections
                // Remember that this function traces only a single ray.
                //   rpl has been set to zero during initialisation.
                float fetchX = (beam_source.x + alpha_mid*diff.x - start.x) / voxel.x;
                float fetchY = (beam_source.y + alpha_mid*diff.y - start.y) / voxel.y;
                float fetchZ = (beam_source.z + alpha_mid*diff.z - start.z) / voxel.z;

                rpl += intersection * tex3D<float>( tex3Dvol, fetchX, fetchY, fetchZ);
            }
            alpha = nextalpha;
        }
        assert(iter<max_iters); // something is wrong
    }
    return rpl;
}

__device__
float cudaSiddon_device(
        float3              this_voxel,
        float*              fluence_map,
        float3              start,
        float3              end,
        float3              voxel,
        float               beamhard_correct,
        float3              beam_source,
        float3              beam_direction,
        float3              beam_isocenter,
        float               beam_sad,
        float2              beamlet_size,
        uint2               fmap_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        int                 dkerns,
        cudaTextureObject_t tex3Ddens,
        cudaTextureObject_t specTex,
        int                 beamletnum
) {
    float terma = 0.f;
    // avoid division by zero (fixes blank center pixel issue)
    /* this_voxel += 1e-9f; */

    // projection vector for this thread
    float3 diff = this_voxel - beam_source;

    /////////////// Check for non-zero fluence /////////////////
    float3 normDir = normalize(beam_direction);
    float3 projection = dot(diff, normDir) * normDir;
    float3 rejection = diff - projection;
    float trad = length(rejection) * fabsf(length(beam_isocenter-beam_source)/length(projection)); // distance from src-axis vector along perpendicular plane
    float3 normRejection = normalize(rejection);

    // calculate conformality based on precomputed fluence map
    float3 proj_fmap = trad * normRejection;
    proj_fmap = rotateBeamAtOriginRHS(proj_fmap, beam_azimuth, beam_zenith, beam_coll);
    float3 fluence_check = proj_fmap;

    // convert coords on fluence map plane to indices
    fluence_check.x = (fluence_check.x/beamlet_size.x) + (fmap_size.x * 0.5f);
    fluence_check.z = (fluence_check.z/beamlet_size.y) + (fmap_size.y * 0.5f);
    if (!(fluence_check.x < 0 || fluence_check.x >= fmap_size.x || fluence_check.z < 0 || fluence_check.z >= fmap_size.y)) {
        int fluence_idx = __float2int_rd(fluence_check.x) + fmap_size.x * (__float2int_rd(fluence_check.z));
        if (!(fluence_idx < 0 || fluence_idx >= fmap_size.x*fmap_size.y || (beamletnum>=0 && fluence_idx != beamletnum))) {
            float fluence = fluence_map[fluence_idx];
            if (!(fluence <= 0.f)) {
                ////////////////////////////////////////////////////////////
                float rpl = cudaRaytrace_device(
                        this_voxel,
                        start,
                        end,
                        voxel,
                        beam_source,
                        tex3Ddens
                        );

                // TODO: add insideness calculation for fraction of voxel in field
                // TODO: alternatively, supersample terma calculation to reduce aliasing
                /* float insideness = 1.f; */

                // use the effective radiological depth of this voxel to find the energy deposited in it by the primary beam
                // dkerns: number of discrete energy levels in beam spectrum
                float tsum=0.f;
                /* float ksum=0.f; */
                for (int e=0; e<dkerns; e++) {
                    float this_fluence = fluence * tex2D<float>( specTex, __int2float_rn(e), 0.0f);
                    float energy = tex2D<float>( specTex, __int2float_rn(e), 1.0f);
                    float mu_en = tex2D<float>( specTex, __int2float_rn(e), 2.0f);
                    float mu = tex2D<float>( specTex, __int2float_rn(e), 3.0f);

                    tsum += this_fluence * energy * mu * exp(-1.f * mu * rpl);
                    /* ksum += this_fluence * energy * mu_en * exp(-1.f * mu * rpl); */
                }

                //Inverse square correction to dose (see Papanikolaou and Mackie 1993)
                // Papanikolaou and Mackie 1993 used inv. sq. correction on dose at dose deposition site
                //    to correct for lack of kernel tilting. Since we handle kernel tilting in the beamlet calculator,
                //    this should instead be used on the TERMA
                float invsq = __powf(beam_sad, 2)/__powf(length(this_voxel - beam_source), 2);

                //calculate T and Kc at zero depth for use in hardening correction
                // beam hardening correction for single poly-energetic kernel
                //see Hoban et al 1994 (PMB)
                terma = tsum * beamhard_correct * invsq;
            }
        }
    }
    return terma;
}

__global__
/* __launch_bounds__(256, 12) // We want to hide global memory latency by decreasin #registers/thread and increasing SM occupancy */
void cudaSiddon(
        float*              terma,
        float*              fluence_map,
        float3              start,
        float3              f_size,
        float3              voxel,
        uint3               calc_bbox_start,
        uint3               calc_bbox_size,
        float               beamhard_correct,
        float3              beam_source,
        float3              beam_direction,
        float3              beam_isocenter,
        float               beam_sad,
        float2              beamlet_size,
        uint2               fmap_size,
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        int                 dkerns,
        cudaTextureObject_t tex3Ddens,
        cudaTextureObject_t specTex,
        int                 beamletnum
) {
    // supersampling setup
    extern __shared__ float s_subterma[];
    s_subterma[s_index(subvoxidx)] = 0.f;
    __syncthreads();


	// launch a thread for each voxel in the cartesian data volume
    int ray_X = threadIdx.y + blockIdx.y * blockDim.y;
    int ray_Y = threadIdx.z + blockIdx.z * blockDim.z;
    int ray_Z = blockIdx.x;

    if (!(ray_X >= calc_bbox_size.x || ray_Y >= calc_bbox_size.y || ray_Z >= calc_bbox_size.z)) {
        // get sub-vox center position
        float3 this_voxel;
        if (subvoxcount > 1) {
            int ss_factor = cbrtf(subvoxcount);
            int sub_X = subvoxidx % ss_factor;
            int sub_Y = (subvoxidx / ss_factor)%ss_factor;
            int sub_Z = subvoxidx / (ss_factor*ss_factor);

            float3 subvoxsize = voxel/(float)ss_factor;
            this_voxel = make_float3(
                    start.x + (calc_bbox_start.x + ray_X-0.5f)*voxel.x + (sub_X+0.5f)*subvoxsize.x,
                    start.y + (calc_bbox_start.y + ray_Y-0.5f)*voxel.y + (sub_Y+0.5f)*subvoxsize.y,
                    start.z + (calc_bbox_start.z + ray_Z-0.5f)*voxel.z + (sub_Z+0.5f)*subvoxsize.z );

        } else {
            // for each voxel, the ray begins at the beam source (point source assumption)
            // the ray ends at the center of current voxel
            this_voxel = make_float3(
                    (start.x + (calc_bbox_start.x + ray_X)*voxel.x),
                    (start.y + (calc_bbox_start.y + ray_Y)*voxel.y),
                    (start.z + (calc_bbox_start.z + ray_Z)*voxel.z) );
        }

        // coords of start and end in RCS
        float3 end;
        end.x = start.x + f_size.x * voxel.x;
        end.y = start.y + f_size.y * voxel.y;
        end.z = start.z + f_size.z * voxel.z;
        s_subterma[s_index(subvoxidx)] = cudaSiddon_device(
                this_voxel,
                fluence_map,
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
    }

// -------------------------
    __syncthreads();
    float avgterma = 0.f;
    if (!(ray_X >= calc_bbox_size.x || ray_Y >= calc_bbox_size.y || ray_Z >= calc_bbox_size.z) && subvoxidx==0) {
        for (int ii=0; ii<subvoxcount; ii++) {
            avgterma += s_subterma[s_index(ii)];
        }
        avgterma /= subvoxcount;
        siddon_write(avgterma);
    }
}
