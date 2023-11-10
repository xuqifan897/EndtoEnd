#ifndef __CUDADOSECALC_H__
#define __CUDADOSECALC_H__

#include "brain_defs.h"
#include "beam.h"
#include "kernel.h"

#include <vector>

namespace old
{
    int radconvolveTexture (
        MONO_KERNELS        *mono,
        CONSTANTS           *constants,
        std::vector<BEAM>&  beams,
        int                 nrays
    );
}

#endif