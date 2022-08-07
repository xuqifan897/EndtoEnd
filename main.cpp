#include <vector>
#include <iostream>
#include <fstream>
#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;


int main(int argc, char** argv)
{
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }
    // deviceProperty();

    // phantom initialization
    phantom Phtm;
    phantom_init_default(Phtm);
    Phtm.to_device();
    Phtm.textureInit();
    // Phtm.textureDecon();

    // beam initialization
    vector<beam> beams;
    beams_init(beams);

    // kernel initialization
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();
    (*kernel).texInit();

    optimize_stationary(beams, Phtm);
}


int main_module_test(int argc, char** argv)
{
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }
    // deviceProperty();

    // the following code block is for testing
    // phantom initialization
    phantom Phtm;
    phantom_init_default(Phtm);
    Phtm.to_device();
    Phtm.textureInit();

    // kernel initialization
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();
    (*kernel).texInit();

    // test_modules(Phtm);
    test_FCBB_water_phantom(Phtm);
}


int main_round2_bcmk(int argc, char** argv)
{
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }
    // deviceProperty();

    // phantom initialization
    phantom Phtm;
    phantom_init_default(Phtm);
    Phtm.to_device();
    Phtm.textureInit();
    // Phtm.textureDecon();

    // beam initialization
    vector<beam> beams;
    beams_init(beams);

    // kernel initialization
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();
    (*kernel).texInit();

    // testConvKernel(FCBB6MeV);
    // test_FCBB_water_phantom(Phtm);
    optimize_stationary(beams, Phtm);
    // optimize_stationary_graph(beams, Phtm);
    // test_dose_sum(beams, Phtm);
    // test_calc_FCBB_PVCS_dose_grad(beams, Phtm);
    // test_FCBB_PVCS_backward(beams, Phtm);
    // test_volume_rendering();
    // test_minus_coordinates_of_texture_memory_out_of_curiosity();
    // test_FCBB_BEV_backward(beams, Phtm);
    // test_element_wise_square();
    // test_fluence_map_update(beams);
}


int main_round1_bcmk(int argc, char** argv)
{
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }
    // deviceProperty();

    // testDepthDose(FCBB6MeV);
    // testConvKernel(FCBB6MeV);

    // phantom Phtm;
    // phantom_init_default(Phtm);
    // Phtm.to_device();
    // Phtm.textureInit();
    // // Phtm.textureDecon();
    // runTest(Phtm);

    // vector<beam> beams;
    // beams_init(beams);
    
    // test_convolve();
    // test_convolveT();

    // test_volume_rendering();
    // test_FCBB_init();

    // test_BEV_dose_forward();
    // test_PVCS_surface();
    // test_PVCS_dose_forward();
}