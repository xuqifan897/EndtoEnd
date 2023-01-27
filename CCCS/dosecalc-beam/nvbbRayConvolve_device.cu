#include "nvbbRayConvolve_device.cuh"
#include <assert.h>

// #define DEBUG_PRINTF

// kernel to trace through density and terma volumes
// resampling along convolution direction
// each gpu thread is one of parallel rays. Each thread steps along BEV sampling ray with steps of
// size delR
#define raytrace_write(arr, val) arr[out_idx] = val;
__global__ void
cudaRaytrace(
        float *bDens,                       // storage for BEV resampled density
        float *bTerm,                       // storage for BEV resampled terma
        float beam_azimuth,
        float beam_zenith,
        float beam_coll,
        float kern_theta,
        float kern_phi,   // kernel rotation angles
        float3 revStart,                    // REV volume limit coords relative to BEV start as origin
        float3 bevStart,                    // BEV volume limit coords in RCS (relative to dicom origin)
        uint3 rev_size,                     // size of resampled REV data volumes
        uint3 bev_size,
        uint3 max_rev_size,
        float3 start,
        float3 voxel,
        float3 rev_voxelsize,
        float3 f_calc_bbox_start,
        float3 f_calc_bbox_end,
        cudaTextureObject_t tex3Ddens,      // XYZ input
        cudaTextureObject_t tex3Dter        // XYZ input
) {
    // thread launched for each voxel in REV space
    // Pay attention: we write the data along the x-axis but this maps to the y-axis in the BEV space
    // some uses of the axes will thus be transposed
    int tX = threadIdx.x + blockIdx.x * blockDim.x;
    int tY = threadIdx.y + blockIdx.y * blockDim.y;
    int tZ = blockIdx.z;
    /* // BEV indices */
    /* int bX = tZ; */
    /* int bY = tX; */
    /* int bZ = tY; */
    // check for out of REV bounds
    if (tX >= rev_size.x || tY >= rev_size.y || tZ >= rev_size.z ) { return; }
    int out_idx = tX + (max_rev_size.x * (tY + (max_rev_size.y * tZ)));

    // Rotate around REV origin to get offset relative to revStart in RCS
    float3 bev_idx = inverseRotateKernelAtOriginRHS(__fmaf_rn(revStart, __frcp_rn(rev_voxelsize), make_float3(tZ, tX, tY)), kern_theta, kern_phi);

    // check for out of BEV bounds
    if (bev_idx.x < 0.f || bev_idx.y < 0.f || bev_idx.z < 0.f
        || bev_idx.x >= __int2float_rn(bev_size.x)
        || bev_idx.y >= __int2float_rn(bev_size.y)
        || bev_idx.z >= __int2float_rn(bev_size.z) )
    {
        raytrace_write(bDens, -3.f);
        /* raytrace_write(bTerm, -3.f); */
        return;
    }

    float3 g_dens_indices =  __fmaf_rn(__fsub_rn(bevStart, start), __frcp_rn(voxel), inverseRotateBeamAtOriginRHS(bev_idx, beam_azimuth, beam_zenith, beam_coll));
    // check if out of calc_bbox bounds
    if (g_dens_indices.x < f_calc_bbox_start.x || g_dens_indices.x >= f_calc_bbox_end.x ||
        g_dens_indices.y < f_calc_bbox_start.y || g_dens_indices.y >= f_calc_bbox_end.y ||
        g_dens_indices.z < f_calc_bbox_start.z || g_dens_indices.z >= f_calc_bbox_end.z )
    {
        raytrace_write(bDens, -2.f);
        /* raytrace_write(bTerm, -2.f); */
        return;
    }

    // get RCS coords relative to scanner origin then convert to RCS indices
    // fetch from textures in RCS
    float valDens = tex3D<float>( tex3Ddens, g_dens_indices.x+0.5f, g_dens_indices.y+0.5f, g_dens_indices.z+0.5f );
    float3 g_term_indices = __fsub_rn(g_dens_indices, f_calc_bbox_start);
    float valTerm = tex3D<float>( tex3Dter, g_term_indices.x+0.5f, g_term_indices.y+0.5f, g_term_indices.z+0.5f );

    raytrace_write(bDens, valDens);
    raytrace_write(bTerm, valTerm);
}

// perform a line convolution in both directions along each ray
// each block is parallel ray and each thread is sample point along ray (spaced by delR)
#define conv_tX threadIdx.x
#define conv_tY blockIdx.y
#define conv_tZ blockIdx.z
__global__ void
RowConvolve(
        float *bDens,                    // REV-trans input
        float *bTerma,                   // REV-trans input
        cudaSurfaceObject_t surfDoseObj, // REV dose output
        float f_kern_wt_idx,                 // kernel theta index (for sampling texKern)
        uint3 rev_size,                  // REV volume dims
        uint3 max_rev_size,
        float delr,
        uint nradii,
        uint ntheta,
        uint nphi,
        float penumbra,
        cudaTextureObject_t kernTex      // kernel weights
) {
	// initialize doseArray to 0 since last convolution
    surf3Dwrite(0.f, surfDoseObj, sizeof(float) * conv_tX, conv_tY, conv_tZ, cudaBoundaryModeTrap);

	// each row of data now corresponds to a single ray, as was sampled in the cudaRaytrace kernel
	// a block of threads is launched for each parallel REV ray with a thread for each sampling point along the ray
    int mem_idx = conv_tX + (max_rev_size.x * (conv_tY + (max_rev_size.y * conv_tZ)));

	// dynamically allocated shared memory
    extern __shared__ float shmBlock[];

    // Move density and terma data into shared memory space on each block
    if (conv_tX < rev_size.x && conv_tY < rev_size.y && conv_tZ < rev_size.z ) {
        shmBlock[conv_tX] = bDens[mem_idx];
        shmBlock[conv_tX + rev_size.x] = bTerma[mem_idx];
    }
    __syncthreads();

    // check for out of REV bounds
    if (conv_tX >= rev_size.x || conv_tY >= rev_size.y || conv_tZ >= rev_size.z ) { return; }

    // outside of bev/bbox intersection
    if (shmBlock[conv_tX] < 0.0f) { return; }

	// initialize variables
    float sum = 0.f;
    float r_eff = 0.f, delr_eff = 0.f;
    float cumval = 0.f, last_cumval = 0.f;
    int rad_start = 0;
    int k1 = conv_tX;
    /* int zeroCount = 0; */
    float last_rad = 0.f;
    float next_rad = KERN_RADII[0];
	// determine convolution direction index, for sampling the poly-energetic dose deposition kernel
    float P = f_kern_wt_idx;

	// accumulate convolution contributions in one direction along the row of data
    // This the "forward" direction, such that an interaction point deposits dose at this thread's assignment
    // voxel by traveling in the forward kernel direction
    do {
        k1--;

		// break loop early if sampling outside active calculation volume
        if (k1<0) { break; }

        //Find effective radiological distance to voxel k1 (heterogeneity correction)
        delr_eff = delr * shmBlock[k1]; // shmBlock[k1] == density[k1]
        r_eff += 0.5f * delr_eff; // move to center of this voxel
        if (r_eff > KERN_RADII[nradii - 1]) { break; }

        //Determine radial kernel index for voxel k1
        while( next_rad < r_eff ) {
            if (rad_start == nradii - 1) { break; }
            last_rad = next_rad;
            next_rad = KERN_RADII[++rad_start];
        }

		// accumulate contribution from voxel k1
        float interp = (r_eff - last_rad)/(next_rad-last_rad);
        if (rad_start == 0) {
            cumval = interp * tex2D<float>( kernTex, 0.5f, P + 0.5f);
        } else {
            cumval = tex2D<float>( kernTex, __int2float_rn(rad_start) - 0.5f + interp, P + 0.5f);
        }
        // CCK dose calculation
        float terma = shmBlock[k1+rev_size.x];
        sum += (cumval - last_cumval) * terma;

        last_cumval = cumval;
        r_eff += 0.5f * delr_eff; // move to end of this voxel
    } while (rad_start < nradii - 1);

	// reset variables
    r_eff = delr_eff = 0.f;
    cumval = last_cumval = 0.f;
    rad_start = 0;
    k1 = conv_tX;
    last_rad  = 0.f;
    next_rad = KERN_RADII[0];

	// perform the same operation in the opposite direction along the row of data
    // This the "reverse" direction, such that an interaction point deposits dose at this thread's assignment
    // voxel by traveling in the backward kernel direction
    do {
        k1++;
        if (k1 >= rev_size.x) { break; }

        //Find effective radiological distance to voxel k1
        delr_eff = delr * shmBlock[k1];
        r_eff += 0.5f * delr_eff;
        if (r_eff > KERN_RADII[nradii - 1]) { break; }

        //Determine radial kernel index for voxel k1
        while( next_rad < r_eff ) {
            if (rad_start == nradii - 1) { break; }
            last_rad = next_rad;
            next_rad = KERN_RADII[++rad_start];
        }

        float interp = (r_eff - last_rad)/(next_rad-last_rad);
        if (rad_start == 0) {
            cumval = interp * tex2D<float>( kernTex, 0.5f, __int2float_rn(ntheta)-0.5f-P);
        } else {
            cumval = tex2D<float>( kernTex, __int2float_rn(rad_start) - 0.5f + interp, __int2float_rn(ntheta)-0.5f-P);
        }
        float terma = shmBlock[k1+rev_size.x];
        sum += (cumval - last_cumval) * terma;

        last_cumval = cumval;
        r_eff += 0.5f * delr_eff;
    } while (rad_start < nradii - 1);

    // normalize by nphi
    sum = __fdiv_rn(sum, nphi);

	// write the summed result out using a surface object
    surf3Dwrite(sum, surfDoseObj, sizeof(float) * conv_tX, conv_tY, conv_tZ, cudaBoundaryModeTrap);
}

// accumulate convolution results back onto Cartesian dose grid
#define rev2xyz_append(val) d_dose[tX + (size.x * (tY + (size.y * tZ)))] += val;
__global__ void
REV2XYZdose(
        float               *d_dose,                           // XYZ dose output
        cudaTextureObject_t texDose,             // BEV dose input
        float               beam_azimuth,
        float               beam_zenith,
        float               beam_coll,
        float               kern_theta,
        float               kern_phi,        // kernel rotation angles
        float3              start,
        float3              voxel,
        float3              rev_voxelsize,
        uint3               size,
        uint3               calc_bbox_start,
        uint3               calc_bbox_end,
        float3              revStart,                       // rev box start indices
        uint3               rev_size,
        float3              bevStart
) {
	// a thread is launched for each voxel in the active calculation volume (calc_bbox:XYZ)
    int tX = calc_bbox_start.x + threadIdx.x + blockIdx.x * blockDim.x;
    int tY = calc_bbox_start.y + threadIdx.y + blockIdx.y * blockDim.y;
    int tZ = calc_bbox_start.z + blockIdx.z;

    if (tX >= calc_bbox_end.x || tY >= calc_bbox_end.y || tZ >= calc_bbox_end.z) { return; }

    // bevStart is coords in RCS with respect to machine origin
    // revStart is coords in REV with respect to bevStart
    float3 g_coords = __fmaf_rn(make_float3(tX,tY,tZ), voxel, start);
    float3 b_indices = rotateBeamAtOriginRHS(__fmul_rn(__fsub_rn(g_coords, bevStart), __frcp_rn(voxel)), beam_azimuth, beam_zenith, beam_coll);
    float3 r_indices = __fmaf_rn(-revStart, __frcp_rn(rev_voxelsize), rotateKernelAtOriginRHS(b_indices, kern_theta, kern_phi));

    // ensure the point is within REV boundaries
    if (r_indices.y >= rev_size.x || r_indices.z >= rev_size.y || r_indices.x >= rev_size.z ||
        r_indices.x < 0.f || r_indices.y < 0.f || r_indices.z < 0.f) { return; }

	// retrieve NVB dose results from texture memory
    // coalesced global mem write
    rev2xyz_append(tex3D<float>(texDose, r_indices.y+0.5f, r_indices.z+0.5f, r_indices.x+0.5f));
}
