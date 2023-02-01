#include "binary_io.h"

#include <cstring>
#include <cstdio>
#include <sys/mman.h>
#include <fcntl.h>
#include <sstream>

#include "server/brain_defs.h"
#include "./io_helpers.h"
#include "./kernel.h"

int load_density( SHM_DATA *data ) {
    char *densfile;
    densfile = new char[1024];
    sprintf(densfile,"%s/density.raw",Paths::Instance()->temp_dir().c_str());

    int filed;
    if ( (filed = open( densfile, O_RDONLY )) == -1) {
        printf("Failed to open density file.\n");
        return false;
    }

    data->density = (float *) mmap( NULL, data->size_data * sizeof(float), PROT_READ, MAP_SHARED, filed, 0);
    if (data->density == NULL) {
        printf("Failed to map density file.\n");
        return false;
    }

    delete [] densfile;
    return 1;
}
int load_data( CONSTANTS *host, SHM_DATA *data, bool verbose ) {
    data->size_data = host->nvoxels();

    if ( load_density( data ) < 0) {
        printf("Failed to load density data.\n");
        return false;
    }

    KERNEL kern;
    std::ostringstream kern_fname;
    kern_fname << Paths::Instance()->cum_kernel_file();
    KERNEL::readFromFile(kern, kern_fname.str(), verbose);
    data->radial_boundary = new float[kern.nradii];
    memcpy(data->radial_boundary, kern.radial_boundary, kern.nradii * sizeof(float));
    data->kernel_array = new float[kern.ntheta * kern.nradii];
    memcpy(data->kernel_array, kern.total_matrix.data(), kern.ntheta * kern.nradii * sizeof(float));
    host->ntheta = kern.ntheta;
    host->nradii = kern.nradii;

    int numangles = host->nphi * kern.ntheta/2;
    load_convolution_theta_angles(&host->conv_theta_deg, numangles);
    load_convolution_phi_angles(&host->conv_phi_deg, numangles);

    return 0;
}
int free_data( CONSTANTS* host, SHM_DATA *data ) {
    munmap( data->density, data->size_data * sizeof(float) );
    munmap( host->conv_theta_deg, host->nphi & host->ntheta/2 * sizeof(float) );
    munmap( host->conv_phi_deg, host->nphi & host->ntheta/2 * sizeof(float) );
    delete [] data->kernel_array;
    delete [] data->radial_boundary;
    return 0;
}

