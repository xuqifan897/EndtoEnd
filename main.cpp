#include <vector>
#include <iostream>
#include <fstream>
#include "args.h"
#include "geom.h"

using namespace E2E;
using namespace std;

int main(int argc, char** argv)
{
    if (args_init(argc, argv))
    {
        cerr << "Argument initialization failure." << endl;
        exit;
    }

    testDepthDose(FCBB6MeV);

    // phantom Phtm;
    // phantom_init_default(Phtm);
    // Phtm.to_device();
    // Phtm.textureInit();
    // Phtm.textureDecon();
    // runTest(Phtm);
}