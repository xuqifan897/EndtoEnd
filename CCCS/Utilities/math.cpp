#include "./math.h"

#include <cmath>

bool closeto(float a, float b, float tolerance) {
    return fabsf(b-a) <= tolerance;
}
