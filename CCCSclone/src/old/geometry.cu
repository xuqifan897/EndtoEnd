#include "geometry.h"
#include "dev_intrinsics.h"

CUDEV_FXN float3 old::inverseRotateBeamAtOriginRHS( const float3& vec, const float& theta, const float& phi, const float& coll ) {
    // invert what was done in forward rotation
    float3 tmp = rotateAroundAxisAtOriginRHS(vec, make_float3(0.f, 1.f, 0.f), -(phi+coll)); // coll rotation + correction
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr);                                   // couch rotation
    return rotateAroundAxisAtOriginRHS(tmp, rotation_axis, theta);                          // gantry rotation
}

CUDEV_FXN float3 old::rotateAroundAxisAtOriginRHS(const float3& p, const float3& r, const float& t) {
    // ASSUMES r IS NORMALIZED ALREADY and center is (0, 0, 0)
    // p - vector to rotate
    // r - rotation axis
    // t - rotation angle
    float sptr, cptr;
    fast_sincosf(t, &sptr, &cptr);
    return make_float3(
            (-r.x*(-r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.x*cptr + (-r.z*p.y + r.y*p.z)*sptr,
            (-r.y*(-r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.y*cptr + (+r.z*p.x - r.x*p.z)*sptr,
            (-r.z*(-r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.z*cptr + (-r.y*p.x + r.x*p.y)*sptr
            );
}