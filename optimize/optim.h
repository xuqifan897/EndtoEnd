#ifndef OPTIM
#define OPTIM

#include <vector>
#include "args.h"
#include "geom.h"

namespace E2E
{
    void optimize_stationary(std::vector<beam>& beams, phantom& Phtm);
}

#endif