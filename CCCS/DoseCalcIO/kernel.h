#ifndef __KERNEL_H__
#define __KERNEL_H__

#include <vector>
#include "H5Cpp.h"

/* KERNELS */
/////////////
#define MAX_KERNELS 30	//max number of monoenergetic kernels to use for creating polyenergetic kernel
#define N_KERNEL_RADII 24  //number of radial increments in mono kernels
#define N_KERNEL_ANGLES 48 //number of angular increments in mono kernels
#define N_KERNEL_CATEGORIES 5

//each kernel file contains 7 entries per voxel, the first five are these categores
enum class KERNEL_CATEGORIES {
    primary,
    first_scatter,
    second_scatter,
    multiple_scatter,
    brem_annih
};

//kernel structure for each monoenergetic kernel and the polyenergetic kernel
class KERNEL {
    public:
        std::string kernel_file;
        int   nradii;
        int   ntheta;
        float radial_boundary[N_KERNEL_RADII+1];
        float angular_boundary[N_KERNEL_ANGLES+1];
        std::vector<float> matrix[N_KERNEL_CATEGORIES];  //kernel values for each category
        std::vector<float> total_matrix;				   //sum of all categories (used for current convolution)

        static int readFromFile(KERNEL& kern, std::string fname, bool verbose=false);
        int writeToFile(std::string fname, bool verbose=false);
        static int _readFromHDF5(KERNEL& beam, H5::Group& h5group);
        int _writeToHDF5(H5::Group& h5group) const;

        void print_info();
        void normalize(bool verbose=false);
};

//macros for accessing kernel values
#define KERNEL(kern_array, i, j) kern_array[i + j*N_KERNEL_RADII]
#define KERNEL_VALUE(kern_ptr,category,rr,tt) \
        (kern_ptr)->matrix[static_cast<int>(category)][(rr) + (tt)*(kern_ptr)->nradii]
#define KERNEL_TOTAL_VALUE(kern_ptr,rr,tt) \
        (kern_ptr)->total_matrix[(rr) + (tt)*(kern_ptr)->nradii]

//info and array of KERNEL structures for monoenergetic kernels
struct MONO_KERNELS {
    int    nkernels;
    std::string spectrum_file;
    float  energy[MAX_KERNELS];
    float  fluence[MAX_KERNELS];
    float  mu[MAX_KERNELS];
    float  mu_en[MAX_KERNELS];
    KERNEL kernel[MAX_KERNELS];
};

// read MONO_KERNELS from original file
int read_kernel(KERNEL *kern);

// Combine mono-energetic kernels into a single poly-energetic kernel and write to file
int make_poly_kernel(MONO_KERNELS *mono, KERNEL *poly, int verbose=false, int debug=false);

// reduce kernel into cumulative-cumulative kernel
int make_cumulative_kernel(KERNEL *kern, int NPHI, int NTHETA, int verbose);

// load binary data volume (.raw filetypes) as memory mapped files
int load_convolution_theta_angles( float **ptr, unsigned int size );
int load_convolution_phi_angles( float **ptr, unsigned int size );

#endif // __KERNEL_H__
