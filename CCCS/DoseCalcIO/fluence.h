#ifndef __FLUENCE_H__
#define __FLUENCE_H__

#include <vector>
#include <string>
#include <helper_cuda.h>
#include <helper_math.h>

// forward declarations
class BEAM;

int load_fluence_map( std::vector<float>& fluence_map, uint2 size, std::string& filename, bool verbose=false);
int write_fluence_map( BEAM& beam, int beam_id, uint2 size, bool verbose=false);

#endif // __FLUENCE_H__
