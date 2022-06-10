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
    phantom Phtm;
    phantom_init_default(Phtm);

    // // for debug purposes
    // cout << Phtm.dimension[0] << ", " << Phtm.dimension[1] << ", " << Phtm.dimension[2] << endl;
    // auto& dim = Phtm.dimension;
    // int size = dim[0] * dim[1];
    // float* slice = (float*)malloc(size*sizeof(float));
    // int z = 100;
    // for (int i=0; i<dim[0]; i++)
    // {
    //     for (int j=0; j<dim[1]; j++)
    //     {
    //         int coord_slice = i * dim[1] + j;
    //         int coord = coord_slice * dim[2] + z;
    //         slice[coord_slice] = Phtm.h_HU[coord];
    //     }
    // }
    // string outputFile{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/slice100.dat"};
    // ofstream output_file(outputFile);
    // output_file.write((char*)slice, size*sizeof(float));
}