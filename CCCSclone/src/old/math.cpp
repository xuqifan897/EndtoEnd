#include "math.h"
#include <cmath>

bool old::closeto(float a, float b, float tolerance) {
    return fabsf(b-a) <= tolerance;
}
