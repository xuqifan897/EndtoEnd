# ifndef argparse_h
# define argparse_h
#include <boost/program_options.hpp>
#include <array>
#include <iostream>
namespace po = boost::program_options;

namespace dev
{
    bool argparse(int argc, char** argv);

    extern po::variables_map vm;
    extern std::vector<float> dicomVolumeStartCoords;

    template<typename T>
    const T& getarg(const std::string& key)
    {   
        if (vm.count(key))
            return vm[key].as<T>();
        std::cerr << "The key " << key << " doesn't exist in the argument list!" << std::endl;
        exit(1);
    }

    const std::vector<float>& getarg(const std::string& key);
}


#endif