#ifndef __KERNEL_H__
#define __KERNEL_H__

#include <vector>
#include <string>
#include "H5Cpp.h"

/* KERNELS */
/////////////
#define MAX_KERNELS 30	//max number of monoenergetic kernels to use for creating polyenergetic kernel
#define N_KERNEL_RADII 24  //number of radial increments in mono kernels
#define N_KERNEL_ANGLES 48 //number of angular increments in mono kernels
#define N_KERNEL_CATEGORIES 5

namespace old
{
    //each kernel file contains 7 entries per voxel, the first five are these categores
    enum class KERNEL_CATEGORIES {
        primary,
        first_scatter,
        second_scatter,
        multiple_scatter,
        brem_annih
    };

    class KERNEL
    {
    public:
        std::string kernel_file;
        int nradii;
        int ntheta;
        float radial_boundary[N_KERNEL_RADII+1];
        float angular_boundary[N_KERNEL_ANGLES+1];
        std::vector<float> matrix[N_KERNEL_CATEGORIES];  // kernel values for each category
        std::vector<float> total_matrix;  // sum of all categories (used for current convolution)

        static int readFromFile(KERNEL& kern, std::string fname, bool verbose=false);
        static int _readFromHDF5(KERNEL& kern, H5::Group& h5group);

        void normalize();
        float& get_value(int k, int i, int j) {return this->matrix[k][i + j*this->nradii];}
        float& get_total_vallue(int i, int j) {return this->total_matrix[i + j*this->nradii];}
    };

    //info and array of KERNEL structures for monoenergetic kernels
    struct MONO_KERNELS
    {
        int    nkernels;
        std::string spectrum_file;
        float  energy[MAX_KERNELS];
        float  fluence[MAX_KERNELS];
        float  mu[MAX_KERNELS];
        float  mu_en[MAX_KERNELS];
        KERNEL kernel[MAX_KERNELS];
    };

    int load_convolution_theta_angles(std::vector<float>& conv_theta_deg, int size);
    int load_convolution_phi_angles(std::vector<float>& conv_phi_deg, int size);
}

#endif