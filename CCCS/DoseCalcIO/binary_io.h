#ifndef __BINARY_IO_H__
#define __BINARY_IO_H__

#include <iostream>
#include <cstring>
#include <cerrno>

#include <helper_math.h>
#include "dosecalc_defs.h"
#include "./paths.h"


// forward declarations
struct SHM_DATA;
struct CONSTANTS;

int load_density( SHM_DATA *data );
int load_data( CONSTANTS *host, SHM_DATA *data, bool verbose=false);
int free_data( CONSTANTS* host, SHM_DATA *data );

// load binary data volume (.raw filetype) through normal methods
template<typename T>
int load_binary_data(T *mat, uint3 count, char *filename, bool verbose=false) {
    char *datafile;
    datafile = new char[1024];
    sprintf(datafile,"%s/%s.raw",Paths::Instance()->temp_dir().c_str(),filename);
    FILE *bin;

    if ( (bin = fopen( datafile, "rb" )) == NULL) {
        printf("load_binary_data() failed with error (%d): %s\n", errno, std::strerror(errno));
        return -1;
    }

    size_t sizer = count.x * count.y * count.z;
    unsigned int entries = fread( (void*)mat, sizeof(float), sizer, bin);

    if ( entries != sizer){
        printf("  Binary file has unexpected size! (%d / %lu)\n", entries, sizer);
        fclose(bin);
        return -2;
    }

    delete [] datafile;
    fclose(bin);
    return 1;
}
template<typename T>
int write_debug_data(T *mat, uint3 count, const char* filename, bool verbose=false) {
    return write_debug_data<T>(mat, make_int3(count), filename, verbose);
}
template<typename T>
int write_binary_data(T *mat, uint3 count, const char* filename, bool verbose=false) {
    return write_binary_data<T>(mat, make_int3(count), filename, verbose);
}
template<typename T>
int write_debug_data(T *mat, int3 count, const char* filename, bool verbose=false) {
    char binfile[1024];
    sprintf(binfile,"%s/%s.raw",Paths::Instance()->temp_dir().c_str(),filename);
    return write_binary_data<T>(mat, count, binfile, verbose);
}
template<typename T>
int write_binary_data(T *mat, int3 count, const char* filename, bool verbose=false) {
    // TODO: convert to c++ ofstream with metadata header/footer
    FILE *binout;

    if (verbose) {
        std::cout << "binary data file name is \"" << filename << "\" with size: ("<<count.x << ", "<<count.y<<", "<<count.z<<")"<<std::endl;
    }

    if ( (binout = fopen(filename,"wb")) == NULL) {
        printf("write_binary_data() failed with error (%d): %s\n", errno, std::strerror(errno));
        return(-1); }

    size_t sizer = count.x * count.y * count.z;

    unsigned int entries = fwrite( (const void*)mat, sizeof(float), sizer, binout);

    if ( entries != sizer){
        printf("  Binary file has unexpected size! (%d / %lu)\n", entries, sizer);
        fclose(binout);
        return -2;
    }

    fclose(binout);
    return(1);
}

#endif // __BINARY_IO_H__










