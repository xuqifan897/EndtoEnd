#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <helper_math.h>
#include "macros.h"

namespace old
{
    CUDEV_FXN float3 inverseRotateBeamAtOriginRHS( const float3& vec, const float& theta, const float& phi, const float& coll );
    CUDEV_FXN float3 rotateAroundAxisAtOriginRHS(const float3& p, const float3& r, const float& t);
}

#endif