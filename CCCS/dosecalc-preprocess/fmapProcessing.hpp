#ifndef __FMAPPROCESSING_H__
#define __FMAPPROCESSING_H__

#include <vector>
#include "helper_cuda.h"
#include "helper_math.h"

int fmap_is_apertureready(const std::vector<float>& fluence_map, uint2 fmap_size);
void fmap_post_apertureready(std::vector<float>& fluence_map, uint2 fmap_size, float fill_val=1.f);

#endif // __FMAPPROCESSING_H__
