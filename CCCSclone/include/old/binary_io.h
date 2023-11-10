#ifndef __BINARY_IO_H__
#define __BINARY_IO_H__

#include "brain_defs.h"

#include <iostream>
#include <cstring>
#include <cerrno>

namespace old
{
    int load_data(CONSTANTS* host, SHM_DATA* data);
    int load_density(SHM_DATA* data);
}

#endif