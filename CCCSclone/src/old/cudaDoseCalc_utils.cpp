#include "cudaDoseCalc.h"
#include "geometry.h"
#include "binary_io.h"
#include <limits>

/* Get RCS coordinates of start/end voxels in REV coordinate system rotated to REV orientation
 * This is critical for converting between REV indices to RCS coords then RCS indices to sample terma/density in RCS
 */
void old::update_extents(float3& currMin, float3& currMax, const float3& thisPt) {
    if (thisPt.x < currMin.x) currMin.x = thisPt.x;
    if (thisPt.y < currMin.y) currMin.y = thisPt.y;
    if (thisPt.z < currMin.z) currMin.z = thisPt.z;
    if (thisPt.x > currMax.x) currMax.x = thisPt.x;
    if (thisPt.y > currMax.y) currMax.y = thisPt.y;
    if (thisPt.z > currMax.z) currMax.z = thisPt.z;
    return;
}

void old::findREV (
    REV_DATA  *rev,
    CONSTANTS *constants,
    BEAM      *beam,
    float3    lim_min,
    float3    lim_max,
    float     theta,         // azimuthal angle of the convolution ray [0->pi]
    float     phi,           // zenithal angle of the convolution ray [0->2pi]
    bool      verbose)
{
    // diff in RCS coords relative to origin at min
    // rotated into BEV space relative to origin at min - used to get corner coords in BEV space of terma bounding box
    float3 b_diff = lim_max-lim_min;
    float3 pbev_min = make_float3(std::numeric_limits<float>::max());
    float3 pbev_max = make_float3(std::numeric_limits<float>::min());

    // for 8 corners of the data volume defined by RCS coordinates of BEV box, find position in REV coord sys
    for (int xoff=0; xoff<2; xoff++) {
        for (int yoff=0; yoff<2; yoff++) {
            for (int zoff=0; zoff<2; zoff++) {
                float3 input = make_float3(b_diff.x*xoff, b_diff.y*yoff, b_diff.z*zoff);
                // apply rotation for the convolution ray and evaluate extents in REV space
                float3 output = rotateKernelAtOriginRHS(input, theta, phi);
                // set beam's eye view extents
                update_extents(pbev_min, pbev_max, output);
    } } }

    // rotate limit coords into RCS space and shift relative to RCS origin
    // remember: Pillar_grid to REV orientation is XYZ -> ZXY and rev size reflects REV orientation
    rev->size = make_uint3(
        static_cast<unsigned int>( ceil(fabsf(pbev_max.y-pbev_min.y) / constants->rev_longspacing) ),
        static_cast<unsigned int>( ceil(fabsf(pbev_max.z-pbev_min.z) / constants->rev_latspacing) ),
        static_cast<unsigned int>( ceil(fabsf(pbev_max.x-pbev_min.x) / constants->rev_latspacing) )
    );

    // shrink this bev
    uint3 oldsize = rev->size;
    float3 adjust = {};
    bool was_shrunk = false;
    if (rev->size.x > constants->max_rev_size.x) {
        adjust.x = constants->rev_longspacing*((int)rev->size.x - (int)constants->max_rev_size.x);
        rev->size.x = constants->max_rev_size.x;
        was_shrunk = true;
    }
    if (rev->size.y > constants->max_rev_size.y) {
        adjust.y = constants->rev_latspacing*((int)rev->size.y - (int)constants->max_rev_size.y);
        rev->size.y = constants->max_rev_size.y;
        was_shrunk = true;
    }
    if (rev->size.z > constants->max_rev_size.z) {
        adjust.z = constants->rev_latspacing*((int)rev->size.z - (int)constants->max_rev_size.z);
        rev->size.z = constants->max_rev_size.z;
        was_shrunk = true;
    }
    if (was_shrunk) { pbev_max -= make_float3(adjust.z, adjust.x, adjust.y); }

    // store limits in REV coordinates relative to rotated lim_min, lim_max
    rev->min_coords = pbev_min;
    rev->max_coords = pbev_max;

    if (verbose) {
        printf(" Theta               :  %5.1f deg\n"            , theta*180/PI );
        printf(" Phi                 :  %5.1f deg\n"            , phi*180/PI );
        printf(" REV->size (PG:YZX)  :    %5d x   %5d x   %5d\n", rev->size.x, rev->size.y, rev->size.z);
        printf(" pbev_min            :  %7.2f x %7.2f x %7.2f\n", rev->min_coords.x, rev->min_coords.y, rev->min_coords.z);
        printf(" pbev_max            :  %7.2f x %7.2f x %7.2f\n", rev->max_coords.x, rev->max_coords.y, rev->max_coords.z);
    }
    if (was_shrunk) {
        throw ArraySizeError(theta, phi, oldsize, rev);
    }
}