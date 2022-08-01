#ifndef OPTIM
#define OPTIM

#include <vector>
#include "args.h"
#include "geom.h"

namespace E2E
{
    void optimize_stationary_graph(std::vector<beam>& beams, phantom& Phtm);
    void optimize_stationary(std::vector<beam>& beams, phantom& Phtm);
    void dose_sum(std::vector<beam>& beams, phantom& Phtm, \
        float** d_result, float** d_sources, cudaStream_t stream=0);
    void reduction(float* d_source, uint size, float* d_out0, \
        float* loss, uint idx, cudaStream_t stream=0);
    void test_dose_sum(std::vector<beam>& beams, phantom& Phtm);
    void element_wise_square(float* d_output, float* d_input, uint size, cudaStream_t stream=0);
    void test_element_wise_square();
    void test_modules(phantom& Phtm);

    // the following functions are for computation verification
    void module_test_convolve(beam& Beam);
    void module_test_host_triliner(phantom& Phtm);
    void module_test_host_linear();
    void module_test_BEV_dose_forward(beam& Beam, phantom& Phtm);
}

#endif