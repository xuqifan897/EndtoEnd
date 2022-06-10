#include <array>
#include <fstream>
#include <iostream>
#include "geom.h"
#include "args.h"
#include <string>
#include <vector>

using namespace E2E;
using namespace std;

phantom::phantom()
{
    this->dimension = array<int, 3>({0, 0, 0});
    this->isocenter = array<float, 3>({0, 0, 0});
    this->voxelSize = 0;
    this->h_HU = nullptr;
    this->h_PTVweight = nullptr;
    this->h_PTVtarget = nullptr;
    this->h_OARweight = nullptr;
    this->h_OARtarget = nullptr;
}

phantom::~phantom()
{
    if (this->h_HU != nullptr)
        free(this->h_HU);
    if (this->h_PTVweight != nullptr)
        free(this->h_PTVweight);
    if (this->h_PTVtarget != nullptr)
        free(this->h_PTVtarget);
    if (this->h_OARweight != nullptr)
        free(this->h_OARweight);
    if (this->h_OARtarget != nullptr)
        free(this->h_OARtarget);
}

int E2E::phantom_init_default(phantom& Phtm)
{
    vector<int> dim = get_args<vector<int>>("phantom-dimension");
    vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
    float voxelsize = get_args<float>("voxel-size");
    string phantom_path = get_args<string>("phantom-path");
    string PTV_weight_path = get_args<string>("PTV-weight-path");
    string PTV_target_path = get_args<string>("PTV-target-path");
    string OAR_weight_path = get_args<string>("OAR-weight-path");
    string OAR_target_path = get_args<string>("OAR-target-path");

    Phtm.dimension = array<int, 3>({dim[0], dim[1], dim[2]});
    Phtm.isocenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});
    Phtm.voxelSize = voxelsize/10;
    int size = dim[0] * dim[1] * dim[2];

    vector<float**> pointers{&(Phtm.h_HU), &(Phtm.h_PTVweight), &(Phtm.h_PTVtarget), \
        &(Phtm.h_OARweight), &(Phtm.h_OARtarget)};
    vector<string> files{phantom_path, PTV_weight_path, PTV_target_path, \
        OAR_weight_path, OAR_target_path};
    for (int i=0; i<5; i++)
    {
        *(pointers[i]) = (float*)malloc(size*sizeof(float));
        ifstream input_file(files[i]);
        if (! input_file.is_open())
        {
            cerr << "Could not open this file: " << files[i] << endl;
            return EXIT_FAILURE;
        }
        input_file.read((char*)(*(pointers[i])), size*sizeof(float));
    }

    // normalize
    for (int i=0; i<size; i++)
        Phtm.h_HU[i] /= 1000;

    return 0;
}