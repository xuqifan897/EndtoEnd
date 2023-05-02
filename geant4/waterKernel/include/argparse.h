#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace wk
{
    int argsInit(int argc, char** argv);
    extern po::variables_map vm;

    template<class T>
    T getArg(std::string key)
    {   
        if (vm.count(key))
            return vm[key].as<T>();
        else
            throw std::runtime_error(std::string("the key \"") 
                + key + std::string("\" not initialized"));
    }
}