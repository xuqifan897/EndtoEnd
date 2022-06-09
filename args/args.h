#ifndef ARGS
#define ARGS

#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <iostream>
namespace po = boost::program_options;

namespace E2E
{
    int args_init(int argc, char** argv);
    template<class T>
    T get_args(std::string key);
    
    class CCCSkernel
    {
    public:
        int num_angles;
        float* angles;
        float* Atheta;
        float* Btheta;
        float* atheta;
        float* btheta;

        bool angles_flag;
        bool Atheta_flag;
        bool Btheta_flag;
        bool atheta_flag;
        bool btheta_flag;

        CCCSkernel(int NA=48);
        ~CCCSkernel();
        CCCSkernel(CCCSkernel& old);
        CCCSkernel(CCCSkernel&& old);
    };

    class FCBBkernel
    {
    public:
        int num_depths;
        float* depths;
        float* doses;
        float A;
        float B;
        float a;
        float b;

        bool depths_flag;
        bool doses_flag;

        FCBBkernel(int ND=400);
        ~FCBBkernel();
        FCBBkernel(FCBBkernel& old);
        FCBBkernel(FCBBkernel&& old);
    };

    int spectrum_init();
    int CCCSkernel_init();
    int FCBBkernel_init();

    extern po::variables_map* args;
    extern std::vector<float>* spectrum_energy;
    extern std::vector<float>* spectrum4MeV;
    extern std::vector<float>* spectrum6MeV;
    extern std::vector<float>* spectrum10MeV;
    extern std::vector<float>* spectrum15MeV;
    extern std::vector<float>* spectrum24MeV;

    extern CCCSkernel* CCCS4MeV;
    extern CCCSkernel* CCCS6MeV;
    extern CCCSkernel* CCCS10MeV;
    extern CCCSkernel* CCCS15MeV;
    extern CCCSkernel* CCCS24MeV;

    extern FCBBkernel* FCBB4MeV;
    extern FCBBkernel* FCBB6MeV;
    extern FCBBkernel* FCBB10MeV;
    extern FCBBkernel* FCBB15MeV;
};

template<class T>
T E2E::get_args(std::string key)
{
    // assert(E2E::args!=nullptr);
    // return (*E2E::args)[key].as<T>();
    if (E2E::args == nullptr)
    {
        std::cerr << "arguments are not initialized!" << std::endl;
        exit;
    }
    try
    {
        return (*E2E::args)[key].as<T>();
    }
    catch(const std::exception& e)
    {
        std::cerr << "the key " << key << " is not initialized" << std::endl;
        exit;
    }
}

#endif