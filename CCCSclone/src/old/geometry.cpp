#include "geometry.h"
#include "dev_intrinsics.cuh"
#include "configure.h"
#include "macros.h"

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

// Rotate vec around center using arbitrary axis
float3 old::rotateAroundAxisRHS(const float3& p, const float3& q, const float3& r, const float& t) {
    // ASSUMES r IS NORMALIZED ALREADY
    // p - vector to rotate
    // q - center point
    // r - rotation axis
    // t - rotation angle
    // non-vectorized version
    //    x,y,z = p.(x,y,z)
    //    a,b,c = q.(x,y,z)
    //    u,v,w = r.(x,y,z)
    //
    //    /* (a*(fast_sq(v)+fast_sq(w)) - u*(b*v + c*w - u*x - v*y - w*z))*(1-fast_cosf(t)) + x*fast_cosf(t) + (-c*v + b*w - w*y + v*z)*fast_sinf(t), */
    //    /* (b*(fast_sq(u)+fast_sq(w)) - v*(a*u + c*w - u*x - v*y - w*z))*(1-fast_cosf(t)) + y*fast_cosf(t) + ( c*u - a*w + w*x - u*z)*fast_sinf(t), */
    //    /* (c*(fast_sq(u)+fast_sq(v)) - w*(a*u + b*v - u*x - v*y - w*z))*(1-fast_cosf(t)) + z*fast_cosf(t) + (-b*u + a*v - v*x + u*y)*fast_sinf(t) */

    float sptr, cptr;
    fast_sincosf(t, &sptr, &cptr);
    return make_float3(
            (q.x*(fast_sq(r.y)+fast_sq(r.z)) - r.x*(q.y*r.y + q.z*r.z - r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.x*cptr + (-q.z*r.y + q.y*r.z - r.z*p.y + r.y*p.z)*sptr,
            (q.y*(fast_sq(r.x)+fast_sq(r.z)) - r.y*(q.x*r.x + q.z*r.z - r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.y*cptr + ( q.z*r.x - q.x*r.z + r.z*p.x - r.x*p.z)*sptr,
            (q.z*(fast_sq(r.x)+fast_sq(r.y)) - r.z*(q.x*r.x + q.y*r.y - r.x*p.x - r.y*p.y - r.z*p.z))*(1-cptr) + p.z*cptr + (-q.y*r.x + q.x*r.y - r.y*p.x + r.x*p.y)*sptr
            );
}
float3 old::rotateAroundAxisAtOriginRHS(const float3& p, const float3& r, const float& t) {
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

// convert RCS coords to BEV coords
float3 old::rotateBeamRHS( const float3& vec, const float3& center, const float& theta, const float& phi, const float& coll ) {
    // first rotate around y-axis by phi+coll then rotate point around z'-axis at center by -theta
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr);                          // couch rotation
    float3 tmp = rotateAroundAxisRHS(vec, center, rotation_axis, -theta);          // gantry rotation
    return rotateAroundAxisRHS(tmp, center, make_float3(0.f, 1.f, 0.f), phi+coll); // coll rotation + correction
}
float3 old::rotateBeamAtOriginRHS( const float3& vec, const float& theta, const float& phi, const float& coll ) {
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr);                          // couch rotation
    float3 tmp = rotateAroundAxisAtOriginRHS(vec, rotation_axis, -theta);          // gantry rotation
    return rotateAroundAxisAtOriginRHS(tmp, make_float3(0.f, 1.f, 0.f), phi+coll); // coll rotation + correction
}
// convert BEV coords to RCS coords
float3 old::inverseRotateBeamRHS( const float3& vec, const float3& center, const float& theta, const float& phi, const float& coll ) {
    // invert what was done in forward rotation
    float3 tmp = rotateAroundAxisRHS(vec, center, make_float3(0.f, 1.f, 0.f), -(phi+coll)); // coll rotation + correction
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr);                                   // couch rotation
    return rotateAroundAxisRHS(tmp, center, rotation_axis, theta);                          // gantry rotation
}
float3 old::inverseRotateBeamAtOriginRHS( const float3& vec, const float& theta, const float& phi, const float& coll ) {
    // invert what was done in forward rotation
    float3 tmp = rotateAroundAxisAtOriginRHS(vec, make_float3(0.f, 1.f, 0.f), -(phi+coll)); // coll rotation + correction
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr);                                   // couch rotation
    return rotateAroundAxisAtOriginRHS(tmp, rotation_axis, theta);                          // gantry rotation
}

// convert bev (beamlet) coords to BEV (beam) coords
float3 old::rotateBeamletRHS( const float3& vec, const float3& center, const float& theta, const float& phi) {
    return inverseRotateBeamletRHS(vec, center, -theta, phi);
}
float3 old::rotateBeamletAtOriginRHS( const float3& vec, const float& theta, const float& phi) {
    return inverseRotateBeamletAtOriginRHS(vec, -theta, phi);
}
// convert BEV (beam) coords to bev (beamlet) coords
float3 old::inverseRotateBeamletRHS( const float3& vec, const float3& center, const float& theta, const float& phi) {
    // first rotate around y-axis by -phi then rotate point around z'-axis at center by theta
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr); // couch rotation
    return rotateAroundAxisRHS(vec, center, rotation_axis, -theta);     // gantry rotation
}
float3 old::inverseRotateBeamletAtOriginRHS( const float3& vec, const float& theta, const float& phi) {
    float sptr, cptr;
    fast_sincosf(-phi, &sptr, &cptr);
    float3 rotation_axis = make_float3(sptr, 0.0f, cptr); // couch rotation
    return rotateAroundAxisAtOriginRHS(vec, rotation_axis, -theta);     // gantry rotation
}

// convert BEV coords to REV coords
float3 old::rotateKernelRHS( const float3& vec, const float3& center, const float& theta, const float& phi ) {
    // similar to beam rotation but by first rotating around y-axis by phi then x'-axis at center by theta
    float sptr, cptr;
    fast_sincosf(phi, &sptr, &cptr);
    float3 untilt = rotateAroundAxisRHS(vec, center, make_float3(cptr, 0.f, -sptr), -theta);
    return rotateAroundAxisRHS(untilt, center, make_float3(0.f, 1.f, 0.f), -phi);
}
float3 old::rotateKernelAtOriginRHS( const float3& vec, const float& theta, const float& phi ) {
    float sptr, cptr;
    fast_sincosf(phi, &sptr, &cptr);
    float3 untilt = rotateAroundAxisAtOriginRHS(vec, make_float3(cptr, 0.f, -sptr), -theta);
    return rotateAroundAxisAtOriginRHS(untilt, make_float3(0.f, 1.f, 0.f), -phi);
}
// convert REV coords to BEV coords
float3 old::inverseRotateKernelRHS( const float3& vec, const float3& center, const float& theta, const float& phi ) {
    // undo what was done by rotateKernelRHS
    float3 roll = rotateAroundAxisRHS(vec, center, make_float3(0.f, 1.f, 0.f), phi);               // kernel roll
    float sptr, cptr;
    fast_sincosf(phi, &sptr, &cptr);
    return rotateAroundAxisRHS(roll, center, make_float3(cptr, 0.f, -sptr), theta);  // kernel tilt
}
float3 old::inverseRotateKernelAtOriginRHS( const float3& vec, const float& theta, const float& phi ) {
    float3 roll = rotateAroundAxisAtOriginRHS(vec, make_float3(0.f, 1.f, 0.f), phi);               // kernel roll
    float sptr, cptr;
    fast_sincosf(phi, &sptr, &cptr);
    return rotateAroundAxisAtOriginRHS(roll, make_float3(cptr, 0.f, -sptr), theta);  // kernel tilt
}

/* each beamlet is defined by the anchor coords corresponding to the points along the beamlet
 * central axis that are nearest and furthest from the beam source, and within the calc_bbox
 * This function returns the start/end anchor GCS coords given the beamlet index
 */
#define FCOMP(x, ii) *(((float*)&x)+ii) // return component of float vector by index
void old::calcBeamletAnchors(
    float3&      start,            // out: start_anchor_coords
    float3&      end,              // out: end_anchor_coords
    float2&      beamletAngles,    // out: beam angles + beamlet divergence angles (.x: azi, .y: .zen)
    float3&      beamletIso,       // out: beamlet isocenter coords
    float3       src,              // beam src coords
    float3       iso,              // beam iso coords
    unsigned int beamlet_idx,      // fluence map linearized index
    float2       beamlet_size,     // [unit: cm]
    uint2        fmap_dims,        // number of beamlets along each axis
    float3       voxelsize,        // [unit: cm]
    float3       density_start,    // coords of start of density volume
    uint3        calc_bbox_start,  // nvoxel offset from density start
    uint3        calc_bbox_size,   // nvoxel size from calc_bbox_start
    float        gantry_angle_rad, // [unit: rad]
    float        couch_angle_rad,  // [unit: rad]
    float        coll_angle_rad,   // [unit: rad]
    int          verbose
) {
    // get beamlet central axis coords on fmap (beamlet-isocenter)
    unsigned int bx = beamlet_idx % fmap_dims.x;
    unsigned int by = beamlet_idx / fmap_dims.x;
    float3 bdiff = make_float3(
            beamlet_size.x * (-0.5f*fmap_dims.x + bx + 0.5f),
            0,
            beamlet_size.y * (-0.5f*fmap_dims.y + by + 0.5f) );
    beamletIso = iso + inverseRotateBeamAtOriginRHS(bdiff, gantry_angle_rad, couch_angle_rad, coll_angle_rad);
    // TODO: CHECK THIS (beamletangles)
    beamletAngles = make_float2(
            acosf(length(iso-src)/length(beamletIso-src)),
            atan2f(bdiff.z, bdiff.x) );

    // get coords of bbox limits
    float3 bbox_begin = density_start + voxelsize * make_float3(calc_bbox_start);
    float3 bbox_end   = bbox_begin + voxelsize * make_float3(calc_bbox_size);

    double distance_buffer = 0.05;

    // evaluate intersection with each of 6 calc_bbox faces
    float3 intersections[6] = {};
    start = make_float3(std::numeric_limits<float>::max());
    end = make_float3(std::numeric_limits<float>::min());
    double min_dist = std::numeric_limits<double>::max();
    double max_dist = std::numeric_limits<double>::min();
    for (int ii=0; ii<3; ii++) { // for x, y, z
        float3 diff_src2b = beamletIso - src;
        double denom = FCOMP(diff_src2b, ii);
        double alpha1 = ( FCOMP(bbox_begin,ii) - FCOMP(src,ii) ) / (denom+1e-12);
        double alpha2 = ( FCOMP(bbox_end,ii)   - FCOMP(src,ii) ) / (denom+1e-12);
        intersections[2*ii  ] = src + alpha1 * diff_src2b;
        intersections[2*ii+1] = src + alpha2 * diff_src2b;
        if (verbose>2) {
            for (int _jj=0; _jj<2; _jj++) {
                float3& _tmp = intersections[2*ii+_jj];
                std::cout << "bbox intersection #"<<(ii*2+_jj) <<": "<<_tmp.x << ", "<<_tmp.y<<", "<<_tmp.z<<")"<< std::endl;
            }
        }
    }

    // check for valid intersection with calc_bbox faces
    for (int ii=0; ii<3; ii++) { // for x, y, x
        for (int jj=0; jj<2; jj++) { // for each intersection (of 2 per axis)
            int idx = 2*ii+jj;
            // given intersection with one dimension, do other two dims occur within bounds of box?
            if (FCOMP(intersections[idx],(ii+1)%3)+distance_buffer >= FCOMP(bbox_begin,(ii+1)%3) && FCOMP(intersections[idx],(ii+1)%3)-distance_buffer <= FCOMP(bbox_end,(ii+1)%3) &&
                FCOMP(intersections[idx],(ii+2)%3)+distance_buffer >= FCOMP(bbox_begin,(ii+2)%3) && FCOMP(intersections[idx],(ii+2)%3)-distance_buffer <= FCOMP(bbox_end,(ii+2)%3) )
            {
                if (verbose>2) {
                    float3& _tmp = intersections[idx];
                    std::cout << "valid intersection found (#"<<idx<<"): "<<_tmp.x << ", "<<_tmp.y<<", "<<_tmp.z<<")" << std::endl;
                }
                float3 _inter = intersections[idx];
                if (length(_inter - src) < min_dist) {
                    start = _inter;
                    min_dist = length(start - src);
                }
                if (length(_inter - src) > max_dist) {
                    end = _inter;
                    max_dist = length(end - src);
                }
            }
        }
    }

    if (!(min_dist <std::numeric_limits<float>::max() && max_dist > std::numeric_limits<float>::min())) {
        std::ostringstream msg;
        msg << "Failed to determine beamlet anchor coordinates for beamlet #" + std::to_string(beamlet_idx);
        if (verbose > 1) {
        msg << std::endl <<
            "DEBUG INFO:" << std::endl <<
            "beamlet angles:  az:"<<beamletAngles.x*180.f/PI <<"; zen:"<<beamletAngles.y*180.f/PI<<" deg)" << std::endl <<
            "bdiff:      ("<<bdiff.x<<","<<bdiff.y<<","<<bdiff.z<<")" << std::endl <<
            "biso:       ("<<beamletIso.x<<","<<beamletIso.y<<","<<beamletIso.z<<")" << std::endl <<
            "start:      ("<<start.x<<","<<start.y<<","<<start.z<<")" << std::endl <<
            "end:        ("<<end.x<<","<<end.y<<","<<end.z<<")" << std::endl <<
            "bbox_begin: ("<<bbox_begin.x<<","<<bbox_begin.y<<","<<bbox_begin.z<<")" << std::endl <<
            "bbox_end:   ("<<bbox_end.x<<","<<bbox_end.y<<","<<bbox_end.z<<")" << std::endl;
        }
        throw std::runtime_error(msg.str());
    }
}