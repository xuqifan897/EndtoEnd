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
    // runTest(Phtm);

    // beam initialization
    vector<beam> beams;
    beams_init(beams);

    // kernel initialization
    FCBBkernel* kernel = FCBB6MeV;
    (*kernel).d_conv_kernel_init();
    (*kernel).texInit();

    optimize_stationary(beams, Phtm);
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
    test_FCBB_water_phantom();
}