#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "args.h"

#define PI 3.141592653589793238

using namespace E2E;
using namespace std;

vector<float>* E2E::spectrum_energy = nullptr;
vector<float>* E2E::spectrum4MeV = nullptr;
vector<float>* E2E::spectrum6MeV = nullptr;
vector<float>* E2E::spectrum10MeV = nullptr;
vector<float>* E2E::spectrum15MeV = nullptr;
vector<float>* E2E::spectrum24MeV = nullptr;

CCCSkernel* E2E::CCCS4MeV = nullptr;
CCCSkernel* E2E::CCCS6MeV = nullptr;
CCCSkernel* E2E::CCCS10MeV = nullptr;
CCCSkernel* E2E::CCCS15MeV = nullptr;
CCCSkernel* E2E::CCCS24MeV = nullptr;

FCBBkernel* E2E::FCBB4MeV = nullptr;
FCBBkernel* E2E::FCBB6MeV = nullptr;
FCBBkernel* E2E::FCBB10MeV = nullptr;
FCBBkernel* E2E::FCBB15MeV = nullptr;

int E2E::spectrum_init()
{
    string spectrum_path = get_args<string>("spectrum-path");
    ifstream input_file(spectrum_path);
    if (!input_file.is_open())
    {
        cerr << "Could not open this file: " << spectrum_path << endl;
        return EXIT_FAILURE;
    }
    string line;
    vector<string> lines;
    while (getline(input_file, line))
        lines.push_back(line);
    input_file.close();
    // for (const auto& i : lines)
    //     cout << i << endl;

    vector<vector<float>**> spectrums {&E2E::spectrum_energy, &E2E::spectrum4MeV, 
        &E2E::spectrum6MeV, &E2E::spectrum10MeV, &E2E::spectrum15MeV, &E2E::spectrum24MeV};
    for (vector<float>** a : spectrums)
        (*a) = new vector<float>(lines.size(), 0);

    for (int i=0; i<lines.size(); i++)
    {
        vector<string> line_parsed;
        stringstream ss(lines[i]);
        string item;
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<spectrums.size(); j++)
        {   
            item = line_parsed[j];
            (**(spectrums[j]))[i] = (item.size()==0) ? 0 : stof(item);
        }
    }

    // for (auto& spectrum : spectrums)
    // {
    //     for (int i=0; i<lines.size(); i++)
    //     {
    //         cout << (**spectrum)[i] << ", ";
    //     }
    //     cout << "\n\n";
    // }

    return 0;
}

int E2E::CCCSkernel_init()
{
    string ATheta_path = get_args<string>("ATheta-path");
    string BTheta_path = get_args<string>("BTheta-path");
    string atheta_path = get_args<string>("atheta-path");
    string btheta_path = get_args<string>("btheta-path");

    ifstream input_file(ATheta_path);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << ATheta_path << endl;
        return EXIT_FAILURE;
    }
    string line;
    vector<string> lines;
    while (getline(input_file, line))
        lines.push_back(line);
    input_file.close();

    int num_angles = 48;
    int num_kernels = 5;
    float angle_base = 1.875 * PI / 180;
    CCCSkernel*** kernels = new CCCSkernel**[num_kernels]\
        {&CCCS4MeV, &CCCS6MeV, &CCCS10MeV, &CCCS15MeV, &CCCS24MeV};
    
    for (int i=0; i<num_kernels; i++)
    {
        *(kernels[i]) = new CCCSkernel(num_angles);
        for (int j=0; j<num_angles; j++)
            (**(kernels[i])).angles[j] = (2*j+1)*angle_base;
    }

    for (int i=0; i<lines.size(); i++)
    {
        stringstream ss(lines[i]);
        string item;
        vector<string> line_parsed;
        while(getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<num_kernels; j++)
            (**(kernels[j])).Atheta[i] = stof(line_parsed[j+1]);
    }

    input_file.open(BTheta_path);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << BTheta_path << endl;
        return EXIT_FAILURE;
    }
    lines.clear();
    while (getline(input_file, line))
        lines.push_back(line);
    input_file.close();

    for (int i=0; i<lines.size(); i++)
    {
        stringstream ss(lines[i]);
        string item;
        vector<string> line_parsed;
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<num_kernels; j++)
            (**(kernels[j])).Btheta[i] = stof(line_parsed[j+1]);
    }

    input_file.open(atheta_path);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << atheta_path << endl;
        return EXIT_FAILURE;
    }
    lines.clear();
    while (getline(input_file, line))
        lines.push_back(line);
    input_file.close();

    for (int i=0; i<lines.size(); i++)
    {
        stringstream ss(lines[i]);
        string item;
        vector<string> line_parsed;
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<num_kernels; j++)
            (**(kernels[j])).atheta[i] = stof(line_parsed[j+1]);
    }

    input_file.open(btheta_path);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << btheta_path << endl;
        return EXIT_FAILURE;
    }
    lines.clear();
    while (getline(input_file, line))
        lines.push_back(line);
    input_file.close();

    for (int i=0; i<lines.size(); i++)
    {
        stringstream ss(lines[i]);
        string item;
        vector<string> line_parsed;
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<num_kernels; j++)
            (**(kernels[j])).btheta[i] = stof(line_parsed[j+1]);
    }

    // // for debug purpose
    // for(int i=0; i<num_kernels; i++)
    // {
    //     CCCSkernel& kernel = **(kernels[i]);
    //     for (int j=0; j<kernel.num_angles; j++)
    //         cout << kernel.angles[j] << ", ";
    //     cout << "\n\n";

    //     for (int j=0; j<kernel.num_angles; j++)
    //         cout << kernel.Atheta[j] << ", ";
    //     cout << "\n\n";

    //     for (int j=0; j<kernel.num_angles; j++)
    //         cout << kernel.Btheta[j] << ", ";
    //     cout << "\n\n";

    //     for (int j=0; j<kernel.num_angles; j++)
    //         cout << kernel.atheta[j] << ", ";
    //     cout << "\n\n";

    //     for (int j=0; j<kernel.num_angles; j++)
    //         cout << kernel.btheta[j] << ", ";
    //     cout << "\n\n";
    // }

    return 0;
}

int E2E::FCBBkernel_init()
{
    string FCBBkernelPath = get_args<string>("pencil-path");
    string depthDosePath = get_args<string>("depthDose-path");

    int num_kernels=4;
    FCBBkernel*** kernels = new FCBBkernel**[num_kernels]{&FCBB4MeV, &FCBB6MeV, &FCBB10MeV, &FCBB15MeV};
    ifstream input_file(depthDosePath);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << depthDosePath << endl;
        return EXIT_FAILURE;
    }
    vector<string> lines;
    string line;
    while(getline(input_file, line))
        lines.push_back(line);
    input_file.close();

    for (int i=0; i<num_kernels; i++)
        (*(kernels[i])) = new FCBBkernel(lines.size());
    
    for (int i=0; i<lines.size(); i++)
    {
        vector<string> line_parsed;
        string item;
        stringstream ss(lines[i]);
        while(getline(ss, item, ','))
            line_parsed.push_back(item);
        for (int j=0; j<num_kernels; j++)
        {
            (**(kernels[j])).depths[i] = stof(line_parsed[0]);
            (**(kernels[j])).doses[i] = stof(line_parsed[j+1]);
        }
    }

    input_file.open(FCBBkernelPath);
    if (! input_file.is_open())
    {
        cerr << "Could not open this file: " << FCBBkernelPath << endl;
        return EXIT_FAILURE;
    }

    lines.clear();
    while(getline(input_file, line))
        lines.push_back(line);
    if (lines.size() != num_kernels)
    {
        cerr << "File error!" << endl;
        return EXIT_FAILURE;
    }
    for (int i=0; i<num_kernels; i++)
    {
        vector<string> line_parsed;
        string item;
        stringstream ss(lines[i]);
        while (getline(ss, item, ','))
            line_parsed.push_back(item);
        FCBBkernel& kernel = (**(kernels[i]));
        kernel.A = stof(line_parsed[0]);
        kernel.B = stof(line_parsed[1]);
        kernel.a = stof(line_parsed[2]);
        kernel.b = stof(line_parsed[3]);
    }

    // // for debug purpose
    // for (int i=0; i<num_kernels; i++)
    // {
    //     FCBBkernel& kernel = (**(kernels[i]));
    //     for (int j=0; j<kernel.num_depths; j++)
    //         cout << kernel.depths[j] << ", " << kernel.doses[j] << endl;
    //     cout << "\n\n";
    // }
    // return 0;

    // // for debug purpose
    // for (int i=0; i<num_kernels; i++)
    // {
    //     FCBBkernel& kernel = (**(kernels[i]));
    //     cout << kernel.A << " " << kernel.B << " " << kernel.a << " " << kernel.b << endl;
    // }
}