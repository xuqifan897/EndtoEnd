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
    void reduction_small(float* d_source, uint size, float* loss, \
        uint idx, cudaStream_t stream=0);
    void test_dose_sum(std::vector<beam>& beams, phantom& Phtm);
    void element_wise_square(float* d_output, float* d_input, uint size, cudaStream_t stream=0);
    void test_element_wise_square();
    void test_modules(phantom& Phtm);
    void fluence_map_init(std::vector<beam>& beams, phantom& Phtm);
    int find_minimum_index(float* pointer, int range);

    // the following functions are for computation verification
    void test_modules_beam_init(beam& Beam, phantom& Phtm);
    void host_BEV_to_PVCS(float PVCS_coords[3], float BEV_coords[3], float theta, float phi);
    void host_PVCS_to_BEV(float BEV_coords[3], float PVCS_coords[3], float theta, float phi);
    void module_test_convolve(beam& Beam);
    float read_phantom(float coords[3], std::array<uint, 3>& phantom_dimension, float* phantom);
    void module_test_host_triliner(phantom& Phtm);
    float read_depth_dose(float* h_depth_dose, float coord, uint size);
    void module_test_host_linear();
    void module_test_BEV_dose_forward(beam& Beam, phantom& Phtm);
    void module_test_BEV_dose_backward(beam& Beam, phantom& Phtm);
    void module_test_PVCS_dose_forward(beam& Beam, phantom& Phtm);
    void write_phantom(float coords[3], std::array<uint, 3>& phantom_dimension, float* phantom, float value);
    void module_test_PVCS_dose_backward(beam& Beam, phantom& Phtm);
    void module_test_dose_sum(std::vector<beam>& beams, phantom& Phtm);
    void module_test_PVCS_dose_grad(phantom& Phtm);
    void module_test_reduction();
    void module_test_fluence_map_update(beam& Beam);
    void module_test_smoothness_calc(beam& Beam, float eta);
    void module_test_small_reduction();
    void module_test_find_minimum_index();
}

#endif