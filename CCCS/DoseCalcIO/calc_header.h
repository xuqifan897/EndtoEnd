#ifndef __CALC_HEADER_H__
#define __CALC_HEADER_H__

#include <helper_cuda.h>
#include <helper_math.h>

#include "./paths.h"

// forward declaration
struct CONSTANTS;

int write_omni_header(
        uint3 count, float3 inc, float3 start,
        uint3 calc_bbox_start, uint3 calc_bbox_size,
        int nphi, int ntheta, int nradii,
        int beam_count, const std::string& beam_spec,
        const std::string& target_structure,
        float rev_latspacing,
        float rev_longspacing,
        float kernel_extent,
        uint ss_factor,
        uint3 max_rev_size,
        float penumbra,
        bool reduce
);
int load_omni_header( CONSTANTS *host, bool verbose );

#endif // __CALC_HEADER_H__
