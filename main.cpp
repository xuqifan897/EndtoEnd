#include <vector>
#include <iostream>
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
    // if (phantom_init_default(Phtm));
    // {
    //     cerr << "Phantom initialization failure." << endl;
    //     exit;
    // }
    phantom_init_default(Phtm);
}