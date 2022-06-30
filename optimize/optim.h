#ifndef OPTIM
#define OPTIM

#include <vector>
#include "args.h"
#include "geom.h"

namespace E2E
{
    void optimize_stationary(std::vector<beam>& beams, phantom& Phtm);
    void dose_sum(std::vector<beam>& beams, phantom& Phtm, \
        float** d_result, cudaStream_t stream=0);
    float reduction(float* d_source, uint size, cudaStream_t stream=0);
    void test_dose_sum(std::vector<beam>& beams, phantom& Phtm);
}

#endif