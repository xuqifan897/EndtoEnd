#ifndef argparse_h
#define argparse_h 1

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace wp
{
    int argsInit(int argc, char** argv);
    extern po::variables_map* vm;
}

#endif