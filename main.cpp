#include <vector>
#include <iostream>
#include "args.h"

using namespace::std;

int main(int argc, char** argv)
{
    E2E::args_init(argc, argv);
    string key = "phantom-path";
    string phantom_path = E2E::get_args<string>(key);
    cout << phantom_path << endl;
}