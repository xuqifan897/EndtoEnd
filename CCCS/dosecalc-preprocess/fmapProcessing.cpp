#include "fmapProcessing.hpp"
#include <cstdio>

struct ActiveSegment {
    int first, last;
};
std::vector<ActiveSegment> get_active_segments(const std::vector<float>& fluence_map, uint2 fmap_size, int yy) {
    auto active = std::vector<ActiveSegment>();
    int begin = -1;
    int end = -1;
    // find active segments
    for (int xx=0; xx<fmap_size.x; ++xx) {
        int idx = fmap_size.x * yy + xx;
        if (fluence_map[idx]!=0) { begin=end=xx; }
        if (begin>=0) {
            if (fluence_map[idx]!=0) { end=xx; }
            if (fluence_map[idx]==0 || xx >= fmap_size.x-1) {
                active.push_back(ActiveSegment{begin, end});
                begin = end = -1;
            }
        }
    }
    return active;
}
int fmap_is_apertureready(const std::vector<float>& fluence_map, uint2 fmap_size) {
    for (int yy=0; yy<fmap_size.y; ++yy) {
        auto active = get_active_segments(fluence_map, fmap_size, yy);
        if (active.size()>1) { return false; }
    }
    return true;
}
void fmap_post_apertureready(std::vector<float>& fluence_map, uint2 fmap_size, float fill_val) {
    /* operate row-by-row, activating any beamlets that violate the "aperture-readyness" of the row.
     * Assignment strategy for activated beamlets is to simply assign '1' */
    for (int yy=0; yy<fmap_size.y; ++yy) {
        auto active = get_active_segments(fluence_map, fmap_size, yy);

        // fill holes
        printf("size: %d\n", active.size());
        for (int ii=0; ii<int(active.size())-1; ++ii) {
            for (int xx=active[ii].last+1; xx<active[ii+1].first; ++xx) {
                fluence_map[fmap_size.x * yy + xx] = fill_val;
            }
        }
    }
}
