#include <array>
#include <fstream>
#include <iostream>
#include "geom.h"
// #include "args.h"
#include <string>
#include <vector>

using namespace E2E;
using namespace std;

phantom::phantom()
{
    this->dimension = array<int, 3>({0, 0, 0});
    this->voxelSize = 0;
    this->isoCenter = array<float, 3>{0, 0, 0};
    this->d_HU = nullptr;
    this->d_PTVweight = nullptr;
    this->d_PTVtarget = nullptr;
    this->d_OARweight = nullptr;
    this->d_OARtarget = nullptr;
}

phantom::~phantom()
{
    if (this->d_HU != nullptr)
        free(this->d_HU);
    if (this->d_PTVweight != nullptr)
        free(this->d_PTVweight);
    if (this->d_PTVtarget != nullptr)
        free(this->d_OARweight);
    if (this->d_OARweight != nullptr)
        free(this->d_OARweight);
    if (this->d_OARtarget != nullptr)
        free(this->d_OARtarget);
}

// int E2E::phantom_init_default(phantom& target)
// {   
//     vector<int> dimension = vector<int> {200, 200, 197};
//     target.voxelSize = 0.2;
//     vector<float> isocenter = vector<float> {205.95534, 211.23352, 162.16011};
//     string HU_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/CT.dat"};
//     string PTVweight_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVweight.dat"};
//     string PTVtarget_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVtarget.dat"};
//     string OARweight_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARweight.dat"};
//     string OARtarget_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARtarget.dat"};

//     // vector<int> dimension = get_args<vector<int>>("phantom-dimension");
//     target.dimension = array<int, 3>({dimension[0], dimension[1], dimension[2]});
//     // target.voxelSize = get_args<float>("voxel-size") / 10; // mm to cm
//     // vector<float> isocenter = get_args<vector<float>>("phantom-isocenter");
//     target.isoCenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});

//     int size = dimension[0] * dimension[1] * dimension[2];
//     target.d_HU = (float*)malloc(size*sizeof(float));
//     target.d_PTVweight = (float*)malloc(size*sizeof(float));
//     target.d_PTVtarget = (float*)malloc(size*sizeof(float));
//     target.d_OARweight = (float*)malloc(size*sizeof(float));
//     target.d_OARtarget = (float*)malloc(size*sizeof(float));

//     // string HU_path = get_args<string>("phantom-path");
//     // string PTVweight_path = get_args<string>("PTV-weight-path");
//     // string PTVtarget_path = get_args<string>("PTV-target-path");
//     // string OARweight_path = get_args<string>("OAR-weight-path");
//     // string OARtarget_path = get_args<string>("OAR-target-path");

//     // float*** pointers = new float**[5]{&(target.d_HU), &(target.d_PTVweight), \
//     //     &(target.d_PTVtarget), &(target.d_OARweight), &(target.d_OARtarget)};target.isoCenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});

//     float*** pointers = (float***)malloc(5*sizeof(float**));
//     pointers[0] = &(target.d_HU);
//     pointers[1] = &(target.d_PTVweight);
//     pointers[2] = &(target.d_PTVtarget);
//     pointers[3] = &(target.d_OARweight);
//     pointers[4] = &(target.d_OARtarget);

//     // vector<string> paths{HU_path, PTVweight_path, PTVtarget_path, OARweight_path, OARtarget_path};
//     // for (int i=0; i<5; i++)
//     // {
//     //     ifstream input_file(paths[i]);
//     //     if (! input_file.is_open())
//     //     {
//     //         cerr << "Could not open this file: " << paths[i] << endl;
//     //         return EXIT_FAILURE;
//     //     }
//     //     input_file.read((char*)(*(pointers[i])), size*sizeof(float));
//     //     input_file.close();
//     // }

//     return 0;
// }

int E2E::phantom_init_default(phantom& target)
{   
    vector<int> dimension = vector<int> {200, 200, 197};
    target.voxelSize = 0.2;
    vector<float> isocenter = vector<float> {205.95534, 211.23352, 162.16011};
    string HU_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/CT.dat"};
    string PTVweight_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVweight.dat"};
    string PTVtarget_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVtarget.dat"};
    string OARweight_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARweight.dat"};
    string OARtarget_path {"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARtarget.dat"};

    target.dimension = array<int, 3>({dimension[0], dimension[1], dimension[2]});
    target.isoCenter = array<float, 3>({isocenter[0]/10, isocenter[1]/10, isocenter[2]/10});

    int size = dimension[0] * dimension[1] * dimension[2];
    target.d_HU = (float*)malloc(size*sizeof(float));
    target.d_PTVweight = (float*)malloc(size*sizeof(float));
    // target.d_PTVtarget = (float*)malloc(size*sizeof(float));
    target.d_OARweight = (float*)malloc(size*sizeof(float));
    target.d_OARtarget = (float*)malloc(size*sizeof(float));

    return 0;
}