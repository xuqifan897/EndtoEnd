#ifndef ARGS
#define ARGS

#include <string>
#include <boost/program_options.hpp>
#include <map>
namespace po = boost::program_options;

namespace E2E
{
    int args_init(int argc, char** argv);
    extern po::variables_map* args;

    template<class T>
    T get_args(std::string key);
};

template<class T>
T E2E::get_args(std::string key)
{
    assert(E2E::args!=nullptr);
    return (*E2E::args)[key].as<T>();
}

#endif