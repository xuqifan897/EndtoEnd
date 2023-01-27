#include "fluence.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "dosecalc_defs.h"
#include "./io_helpers.h"
#include "./beam.h"
#include "./paths.h"

#define FMAP_NAME "fmap-"

using namespace dcio;

int load_fluence_map( std::vector<float>& fluence_map, uint2 size, std::string& filename, bool verbose)
{
    std::ostringstream beamfile;
    beamfile << Paths::Instance()->temp_dir() << "/fluence_maps/" << filename;
    if (verbose) { printf("binary data file name is %s\n", beamfile.str().c_str()); }
    fflush(stdout);

    FILE *bin;
    if ( (bin = fopen(beamfile.str().c_str(),"rb")) == NULL)
    {
        printf("Cannot open binary file\n");
        return(-1);
    }

    size_t sizer = size.x * size.y;

    int entries = fread( (void*)fluence_map.data(), sizeof(float), sizer, bin);
    fclose(bin);

    if ( entries != sizer) {
        printf("  Binary file has unexpected size! (%d / %d)\n", entries, (int)sizer);
        return -2;
    }

    return(1);
}

int write_fluence_map( BEAM& beam, int beam_id, uint2 size, bool verbose)
{
    std::ostringstream fname;
    fname << FMAP_NAME << std::setfill('0') << std::setw(6) << beam_id << std::setw(0) << ".raw";
    beam.fmap_fname = fname.str();
    // sprintf(beam.fmap_fname,"%s%06d.raw",FMAP_NAME, beam_id);

    char dirpath[1024];
    sprintf(dirpath,"%s/fluence_maps",Paths::Instance()->temp_dir().c_str());
    create_directory(dirpath);
    char beamfile[1024];
    sprintf(beamfile,"%s/%s",dirpath,beam.fmap_fname.c_str());
    if (verbose) { printf("binary data file name is %s\n",beamfile); }
    fflush(stdout);

    FILE *binout;
    if ( (binout = fopen(beamfile,"wb")) == NULL) {
        std::cout << "Cannot open binary file" << std::endl;
        return(-1);
    }

    size_t sizer = size.x * size.y;

    int entries = fwrite( (const void*)beam.fluence_map.data(), sizeof(float), sizer, binout);
    fclose(binout);

    if ( entries != (int)sizer) {
        std::cout << "  Binary file has unexpected size! (" << entries << " / " << sizer << ")" << std::endl;
        return -2;
    }

    return(1);
}
