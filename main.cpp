#include <vector>
#include <iostream>
#include "args.h"

using namespace::std;

int main(int argc, char** argv)
{
    if (E2E::args_init(argc, argv)) exit;
    // bool flag = E2E::args_init(argc, argv);
    // if (flag) exit;
    // string key = "phantom-path";
    // string phantom_path = E2E::get_args<string>(key);
    // cout << phantom_path << endl;
}