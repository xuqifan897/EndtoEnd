#include "args.h"
#include "geom.h"
#include "optim.h"

using namespace E2E;
using namespace std;

// int main_module_test(int argc, char** argv)
int main(int argc, char** argv)
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