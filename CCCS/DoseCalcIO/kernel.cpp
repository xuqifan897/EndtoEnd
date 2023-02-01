#include "kernel.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>
#include <vector>

#include "io_helpers.h"
#include "./paths.h"

#include "dosecalc_defs.h"
#include "server/brain_defs.h"
#include "DoseCalcIO/binary_io.h"

// fillers for these entries in kernel structure
#define UNCERT 0.0
#define MEAN_RADIUS 0.0
#define MEAN_ANGLE 0.0

using namespace dcio;

int read_kernel(KERNEL *kern)
{
    int i, j;
    float uncert;
    float mean_radius, mean_angle;

    char kernel_path[350];

    KERNEL_CATEGORIES category;

    FILE *kernfile;
    FILE *radfile;
    FILE *polarfile;


    Paths* paths = Paths::Instance();
    sprintf(kernel_path,"%s/%s",paths->kernel_dir().c_str(),kern->kernel_file.c_str());
    //printf("\nKernel file path is %s\n",kernel_path);

    if ( (kernfile = fopen(kernel_path,"r")) == NULL){
        printf("Cannot open infile %s\n",kern->kernel_file.c_str());
        exit(1);}

    if ( (radfile = fopen(paths->radius_file().c_str(), "r")) == NULL){
        printf("Cannot open infile %s\n",paths->radius_file().c_str());
        exit(1);}

    if ( (polarfile = fopen(paths->polar_file().c_str(),"r")) == NULL){
        printf("Cannot open infile %s\n",paths->polar_file().c_str());
        exit(1);}


    kern->nradii = N_KERNEL_RADII;
    kern->ntheta = N_KERNEL_ANGLES;

    //read radial boundaries
    for (i=0;i<kern->nradii;i++) {
        float temp;
        fscanf(radfile,"%f", &temp);
        kern->radial_boundary[i] = temp;
    }
    fclose(radfile);

    // read mean angles
    float temp_polar;
    for (i=0;i<kern->ntheta;i++) {
        fscanf(polarfile,"%f",&temp_polar);
        // kern->angular_boundary[i] = 180.f * temp_polar/PI;
        kern->angular_boundary[i] = float(i+1)*180.f/(float)kern->ntheta;
    }
    fclose(polarfile);


    //allocate space for each category	and total kernel
    for (i=0;i<N_KERNEL_CATEGORIES;i++)
        kern->matrix[i] = std::vector<float>(kern->ntheta*kern->nradii);
    kern->total_matrix = std::vector<float>(kern->ntheta*kern->nradii);

    float normkern=0;
    //for each kernel element, read values for all categories
    //at each category increment total kernel value for voxel i,j
    for (j=0;j<kern->ntheta;j++)
        for (i=0;i<kern->nradii;i++)
        {
            KERNEL_TOTAL_VALUE(kern,i,j) = 0.0;

            category = KERNEL_CATEGORIES::primary;
            fscanf(kernfile,"%f %f",&(KERNEL_VALUE(kern,category,i,j)),&uncert);
            KERNEL_TOTAL_VALUE(kern,i,j) += KERNEL_VALUE(kern,category,i,j);

            category = KERNEL_CATEGORIES::first_scatter;
            fscanf(kernfile,"%f %f",&(KERNEL_VALUE(kern,category,i,j)),&uncert);
            KERNEL_TOTAL_VALUE(kern,i,j) += KERNEL_VALUE(kern,category,i,j);

            category = KERNEL_CATEGORIES::second_scatter;
            fscanf(kernfile,"%f %f",&(KERNEL_VALUE(kern,category,i,j)),&uncert);
            KERNEL_TOTAL_VALUE(kern,i,j) += KERNEL_VALUE(kern,category,i,j);

            category = KERNEL_CATEGORIES::multiple_scatter;
            fscanf(kernfile,"%f %f",&(KERNEL_VALUE(kern,category,i,j)),&uncert);
            KERNEL_TOTAL_VALUE(kern,i,j) += KERNEL_VALUE(kern,category,i,j);

            category = KERNEL_CATEGORIES::brem_annih;
            fscanf(kernfile,"%f %f",&(KERNEL_VALUE(kern,category,i,j)),&uncert);
            KERNEL_TOTAL_VALUE(kern,i,j) += KERNEL_VALUE(kern,category,i,j);

            normkern += KERNEL_TOTAL_VALUE(kern,i,j);

            //mean radius and angle values - not validated in kernel files (unused)
            fscanf(kernfile,"%f %f",&mean_radius,&uncert);
            fscanf(kernfile,"%f %f",&mean_angle,&uncert);
        }
    fclose(kernfile);

    return(true);
}
int make_poly_kernel(MONO_KERNELS *mono, KERNEL *poly, int verbose, int debug) {
    KERNEL_CATEGORIES category;
    int i, j, e;
    float sum;


    poly->nradii = N_KERNEL_RADII;
    poly->ntheta = N_KERNEL_ANGLES;

    //copy radial boundaries from first mono kernel
    for (i=0;i<poly->nradii;i++)
        poly->radial_boundary[i] = mono->kernel[0].radial_boundary[i];

    //copy  mean angles from first mono kernel
    for (i=0;i<poly->ntheta;i++)
        poly->angular_boundary[i] = mono->kernel[0].angular_boundary[i];

    for (i=0;i<N_KERNEL_CATEGORIES;i++)
        poly->matrix[i] = std::vector<float>(poly->ntheta*poly->nradii);
    poly->total_matrix = std::vector<float>(poly->ntheta*poly->nradii);

    // average value: V_k(Radial, Theta) for all "k" kernel energies
    // see "Effective Mean Kernel Method" (Liu, 1997)
    for (j=0;j<poly->ntheta;j++) {
        for (i=0;i<poly->nradii;i++) {
            KERNEL_TOTAL_VALUE(poly,i,j) = 0.0;

            //weight of each mono kernel value in sum is fluence*energy*mu
            category = KERNEL_CATEGORIES::primary;
            KERNEL_VALUE(poly,category,i,j) = 0.0;
            sum = 0.0;
            for (e=0;e<mono->nkernels;e++)
            {
                KERNEL_VALUE(poly,category,i,j) += mono->fluence[e]*mono->energy[e]*mono->mu[e]
                    * KERNEL_VALUE(&(mono->kernel[e]),category,i,j);
                sum += mono->fluence[e]*mono->energy[e]*mono->mu[e];
            }
            KERNEL_VALUE(poly,category,i,j) /= sum;
            KERNEL_TOTAL_VALUE(poly,i,j) += KERNEL_VALUE(poly,category,i,j);

            category = KERNEL_CATEGORIES::first_scatter;
            KERNEL_VALUE(poly,category,i,j) = 0.0;
            sum = 0.0;
            for (e=0;e<mono->nkernels;e++)
            {
                KERNEL_VALUE(poly,category,i,j) += mono->fluence[e]*mono->energy[e]*mono->mu[e]
                    * KERNEL_VALUE(&(mono->kernel[e]),category,i,j);
                sum += mono->fluence[e]*mono->energy[e]*mono->mu[e];
            }
            KERNEL_VALUE(poly,category,i,j) /= sum;
            KERNEL_TOTAL_VALUE(poly,i,j) += KERNEL_VALUE(poly,category,i,j);

            category = KERNEL_CATEGORIES::second_scatter;
            KERNEL_VALUE(poly,category,i,j) = 0.0;
            sum = 0.0;
            for (e=0;e<mono->nkernels;e++)
            {
                KERNEL_VALUE(poly,category,i,j) += mono->fluence[e]*mono->energy[e]*mono->mu[e]
                    * KERNEL_VALUE(&(mono->kernel[e]),category,i,j);
                sum += mono->fluence[e]*mono->energy[e]*mono->mu[e];
            }
            KERNEL_VALUE(poly,category,i,j) /= sum;
            KERNEL_TOTAL_VALUE(poly,i,j) += KERNEL_VALUE(poly,category,i,j);

            category = KERNEL_CATEGORIES::multiple_scatter;
            KERNEL_VALUE(poly,category,i,j) = 0.0;
            sum = 0.0;
            for (e=0;e<mono->nkernels;e++)
            {
                KERNEL_VALUE(poly,category,i,j) += mono->fluence[e]*mono->energy[e]*mono->mu[e]
                    * KERNEL_VALUE(&(mono->kernel[e]),category,i,j);
                sum += mono->fluence[e]*mono->energy[e]*mono->mu[e];
            }
            KERNEL_VALUE(poly,category,i,j) /= sum;
            KERNEL_TOTAL_VALUE(poly,i,j) += KERNEL_VALUE(poly,category,i,j);

            category = KERNEL_CATEGORIES::brem_annih;
            KERNEL_VALUE(poly,category,i,j) = 0.0;
            sum = 0.0;
            for (e=0;e<mono->nkernels;e++)
            {
                KERNEL_VALUE(poly,category,i,j) += mono->fluence[e]*mono->energy[e]*mono->mu[e]
                    * KERNEL_VALUE(&(mono->kernel[e]),category,i,j);
                sum += mono->fluence[e]*mono->energy[e]*mono->mu[e];
            }
            KERNEL_VALUE(poly,category,i,j) /= sum;
            KERNEL_TOTAL_VALUE(poly,i,j) += KERNEL_VALUE(poly,category,i,j);
        }
    }


    if (debug) {
        FILE *poly_file;

        //poly kernel file is named after spectrum file
        std::ostringstream sstream;
        sstream << Paths::Instance()->temp_dir() << "/" << mono->spectrum_file << ".kern";
        poly->kernel_file = sstream.str();

        if ( (poly_file = fopen(poly->kernel_file.c_str(),"w")) == NULL){
            printf("Cannot open polyenergetic kernel file for writing!\n");
            return(-1);}

        //write out polyenergetic kernel file for debug purposes (not used downstream, since CCK is instead)
        for (j=0;j<poly->ntheta;j++)
            for (i=0;i<poly->nradii;i++)
            {
                category = KERNEL_CATEGORIES::primary;
                fprintf(poly_file," %10.6e   %10.6e\n",KERNEL_VALUE(poly,category,i,j),UNCERT);

                category = KERNEL_CATEGORIES::first_scatter;
                fprintf(poly_file," %10.6e   %10.6e\n",KERNEL_VALUE(poly,category,i,j),UNCERT);

                category = KERNEL_CATEGORIES::second_scatter;
                fprintf(poly_file," %10.6e   %10.6e\n",KERNEL_VALUE(poly,category,i,j),UNCERT);

                category = KERNEL_CATEGORIES::multiple_scatter;
                fprintf(poly_file," %10.6e   %10.6e\n",KERNEL_VALUE(poly,category,i,j),UNCERT);

                category = KERNEL_CATEGORIES::brem_annih;
                fprintf(poly_file," %10.6e   %10.6e\n",KERNEL_VALUE(poly,category,i,j),UNCERT);

                // fprintf(poly_file," %10.6e   %10.6e\n",MEAN_ANGLE,UNCERT);
                // fprintf(poly_file," %10.6e   %10.6e\n",MEAN_RADIUS,UNCERT);
                fprintf(poly_file," %10.6e   %10.6e\n",poly->angular_boundary[j],UNCERT);
                fprintf(poly_file," %10.6e   %10.6e\n",poly->radial_boundary[i],UNCERT);
            }
        fclose(poly_file);
        if (verbose) { printf("Poly kernel file written to %s\n",poly->kernel_file.c_str()); }
    }

    // option to write kernel data out to command line
    if (verbose) { poly->print_info(); }

    return(1);
}

int make_cumulative_kernel(KERNEL *kern, int NPHI, int NTHETA, int verbose)
{
    KERNEL cumkern;

    // init
    cumkern.nradii = kern->nradii;
    cumkern.ntheta = NTHETA;	//must divide evenly into N_KERNEL_ANGLES
    for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) { cumkern.matrix[cc] = std::vector<float>(cumkern.ntheta * cumkern.nradii); }
    cumkern.total_matrix = std::vector<float>(cumkern.ntheta * cumkern.nradii);


    // set up new boundaries to perform interpolation based on requested ntheta
    // original code could only split if angular boundaries aligned for smaller ntheta. This performs linear interpolation
    // to split for any ntheta
    for (int tt=0; tt<cumkern.ntheta+1; tt++) {
        cumkern.angular_boundary[tt] = (180.f/cumkern.ntheta)*(tt+1);
    }
    memcpy(cumkern.radial_boundary, kern->radial_boundary, kern->nradii*sizeof(float));

    if (verbose) {
        printf("\nCumulative Kernel Structure:\n");
        printf("  angular boundaries [deg]: {");
        for (int tt=0; tt<cumkern.ntheta; tt++) {
            printf("%g", cumkern.angular_boundary[tt]);
            if (tt < cumkern.ntheta-1) { printf(", "); }
        }
        printf("}\n");
        printf("  radial boundaries [cm]:   {");
        for (int rr=0; rr<cumkern.nradii; rr++) {
            printf("%g", cumkern.radial_boundary[rr]);
            if (rr < cumkern.nradii-1) { printf(", "); }
        }
        printf("}\n\n");
    }

/*
 *     // calculate original boundaries from original mean angles
 *     // Original code just generates angular boundaries as i+1*180/ntheta, then takes midpoint as mean angle during convolution
 *     // We instead store the mean angles in "polar.dat" file and recover the angular boundaries here.
 *     auto all_angular_boundaries = std::vector<float>(kern->ntheta+1);
 *     for (int tt=0; tt<kern->ntheta-1; tt++) {
 *         // all_angular_boundaries[tt] = 0.5*(kern->angular_boundary[tt+1]+kern->angular_boundary[tt]);
 *         all_angular_boundaries[tt] = kern->angular_boundary[tt]; //FIXME
 *     }
 *     // ugly but necessary - mono-kernels don't extend mean angles beyond 180 degs
 *     float dangle = all_angular_boundaries[kern->ntheta-2] - all_angular_boundaries[kern->ntheta-3];
 *     all_angular_boundaries[kern->ntheta-1] = all_angular_boundaries[kern->ntheta-2] + dangle;
 *     all_angular_boundaries[kern->ntheta] = all_angular_boundaries[kern->ntheta-1] + dangle;

 *     auto cum_mean_angles = std::vector<float>(cumkern.ntheta+1);
 *     for (int tt=0; tt<cumkern.ntheta+1; tt++) {
 *         cum_mean_angles[tt] = 0.5f*(cumkern.angular_boundary[tt] + (tt>0 ? cumkern.angular_boundary[tt-1] : 0));
 *     }
 *
 *     if (verbose) {
 *         printf("\nMono Kernel Structure:\n");
 *         printf("  boundary angles [deg]:      {");
 *         for (int tt=0; tt<kern->ntheta+1; tt++) {
 *             printf("%g", all_angular_boundaries[tt]);
 *             if (tt < kern->ntheta+1) { printf(", "); }
 *         }
 *         printf("}\n");
           printf("  angular bounds [deg]:   {");
           for (int tt=0; tt<cumkern.ntheta; tt++) {
               printf("%g", cum_mean_angles[tt]);
               if (tt < cumkern.ntheta-1) { printf(", "); }
           }
           printf("}\n");
 *     }
 *
 *     // Combine collapsed cones to reduce number of angular kernel coefficients
 *     int jj1 = 0; // index of first bin in sum
 *     int jj2 = 0; // index of last bin in sum
 *     float A = 1; // fraction of first bin in sum
 *     float B = 0; // fraction of last bin in sum
 *     for (int tt=0; tt<cumkern.ntheta; tt++) {
 *         // use linear interpolation for new boundaries that fall between existing boundaries
 *         // find index of next boundary greater/equal to new boundary "j"
 *         while (all_angular_boundaries[jj2]<cum_angular_boundaries[tt]) { jj2++; }
 *         B = (cum_angular_boundaries[tt]-all_angular_boundaries[jj2-1])/(all_angular_boundaries[jj2]-all_angular_boundaries[jj2-1]);
 *
 *         for (int rr=0; rr<cumkern.nradii; rr++) {
 *             // perform left end (partial) sum
 *             for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
 *                 KERNEL_VALUE(&cumkern, cc, rr, tt) += A * KERNEL_VALUE(kern, cc, rr, jj1);
 *             }
 *             KERNEL_TOTAL_VALUE(&cumkern, rr, tt) += A * KERNEL_TOTAL_VALUE(kern, rr, jj1);
 *
 *             // perform central sum
 *             for (int jjiter=jj1+1; jjiter<jj2; jjiter++) {
 *                 for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
 *                     KERNEL_VALUE(&cumkern, cc, rr, tt) += KERNEL_VALUE(kern, cc, rr, jjiter);
 *                 }
 *                 KERNEL_TOTAL_VALUE(&cumkern, rr, tt) += KERNEL_TOTAL_VALUE(kern, rr, jjiter);
 *             }
 *
 *             // perform right end (partial) sum
 *             for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
 *                 KERNEL_VALUE(&cumkern, cc, rr, tt) += B * KERNEL_VALUE(kern, cc, rr, jj2);
 *             }
 *             KERNEL_TOTAL_VALUE(&cumkern, rr, tt) += B * KERNEL_TOTAL_VALUE(kern, rr, jj2);
 *         }
 *         if (verbose>1) {
 *             printf("A: %g; B: %g; j: %d; jj1: %d; jj2: %d; Tj: %g; Tjj1: %g; Tjj2: %g\n",
 *                     A, B, tt, jj1, jj2, cum_angular_boundaries[tt], all_angular_boundaries[jj1], all_angular_boundaries[jj2]);
 *         }
 *
 *         // update for next new boundary
 *         A = (1.f-B);
 *         jj1 = jj2;
 *     }
 */

    // aggregate polar samples evenly (ntheta must factor evenly into 48: (1, 2, 4, 6, 8, 12, 24, 48))
    for (int tt=0; tt<cumkern.ntheta; tt++) {
        int mark = tt * kern->ntheta/cumkern.ntheta;
        for (int kk=0; kk<kern->ntheta/cumkern.ntheta; kk++) {
            for (int rr=0; rr<cumkern.nradii; rr++) {
                for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
                    KERNEL_VALUE(&cumkern,cc,rr,tt) += KERNEL_VALUE(kern,cc,rr,mark+kk);
                }
                KERNEL_TOTAL_VALUE(&cumkern,rr,tt) += KERNEL_TOTAL_VALUE(kern,rr,mark+kk);
            }
        }
    }


    //Make cumulative kernel (along radial axis)
    //this is what is used for the dose calculation
    for (int tt=0; tt<cumkern.ntheta; tt++) {
        for (int rr=0; rr<cumkern.nradii; rr++) {
            for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
                if (rr > 0) { KERNEL_VALUE(&cumkern,cc,rr,tt) += KERNEL_VALUE(&cumkern,cc,rr-1,tt); }
            }
            if (rr > 0) { KERNEL_TOTAL_VALUE(&cumkern,rr,tt) += KERNEL_TOTAL_VALUE(&cumkern,rr-1,tt); }
        }
    }

    // determine optimal theta/phi pairings
    int numangles = cumkern.ntheta/2 * NPHI;
    auto conv_theta_deg = std::vector<float>(numangles);
    auto conv_phi_deg = std::vector<float>(numangles);
    for (int pp=0; pp<NPHI; pp++) {
        for (int tt=0; tt<cumkern.ntheta/2; tt++) {
            int aa = pp*cumkern.ntheta/2 + tt;
            conv_theta_deg[aa] = 0.5f*(cumkern.angular_boundary[tt] + (tt>0 ? cumkern.angular_boundary[tt-1] : 0));
            conv_phi_deg[aa] = (pp+0.5)*360. / NPHI;
        }
    }
    if (verbose) {
        printf("Convolution Angles [deg]:\n");
        for (int pp=0; pp<NPHI; pp++) {
            printf("  phi=%5g :: theta={", conv_phi_deg[pp*cumkern.ntheta/2] );
            for (int tt=0; tt<cumkern.ntheta/2; tt++) {
                printf("%8.3f", conv_theta_deg[tt + pp*cumkern.ntheta/2]);
                if (tt < cumkern.ntheta/2 - 1) { printf(", "); }
            }
            printf("}\n");
        }
        printf("\n");
    }

    //TODO replace with h5 kernel format
    Paths* paths = Paths::Instance();
    write_debug_data<float>( conv_theta_deg.data(), make_uint3( numangles, 1 , 1 ), paths->conv_theta_file().c_str());
    write_debug_data<float>( conv_phi_deg.data(), make_uint3( numangles, 1 , 1 ), paths->conv_phi_file().c_str());
    std::ostringstream fname;
    fname << paths->cum_kernel_file();
    cumkern.writeToFile(fname.str(), verbose);

    return(true);
}

int _load_conv_angles(float **ptr, unsigned int size, char* fname ) {

    FILE* filed;
    if ( (filed = fopen( fname, "rb" )) == NULL)
    {
        printf("Failed to open convolution angles file.\n");
        return false;
    }

    *ptr = new float[size];
    fread(*ptr, sizeof(float), size, filed);
    fclose(filed);

    // *ptr = (float *) mmap( NULL, size * sizeof(float), PROT_READ, MAP_SHARED, filed, 0);
    return true;
}
int load_convolution_theta_angles( float **ptr, unsigned int size ) {
    char *f;
    f = new char[1024];
    sprintf(f, "%s/%s.raw",Paths::Instance()->temp_dir().c_str(), Paths::Instance()->conv_theta_file().c_str());
    bool result =_load_conv_angles(ptr, size, f);
    delete [] f;
    return result;
}
int load_convolution_phi_angles( float **ptr, unsigned int size ) {
    char *f;
    f = new char[1024];
    sprintf(f, "%s/%s.raw",Paths::Instance()->temp_dir().c_str(),Paths::Instance()->conv_phi_file().c_str());
    bool result =_load_conv_angles(ptr, size, f);
    delete [] f;
    return result;
}

void KERNEL::print_info() {
    float integral[N_KERNEL_CATEGORIES] {0.f};
    float total_integral = 0;

    //small kernel used in calc_dose has no name
    if (kernel_file.length() > 1) { printf("\n\nKernel file: %s\n",kernel_file.c_str()); }
    else { printf("\n\nKernel file: Not named\n"); }

    //at each i,j increment integral for each category
    KERNEL_CATEGORIES category;
    for (int j=0;j<ntheta;j++) {
        for (int i=0;i<nradii;i++) {
            category = KERNEL_CATEGORIES::primary;
            integral[static_cast<int>(category)] += KERNEL_VALUE(this,category,i,j);

            category = KERNEL_CATEGORIES::first_scatter;
            integral[static_cast<int>(category)] += KERNEL_VALUE(this,category,i,j);

            category = KERNEL_CATEGORIES::second_scatter;
            integral[static_cast<int>(category)] += KERNEL_VALUE(this,category,i,j);

            category = KERNEL_CATEGORIES::multiple_scatter;
            integral[static_cast<int>(category)] += KERNEL_VALUE(this,category,i,j);

            category = KERNEL_CATEGORIES::brem_annih;
            integral[static_cast<int>(category)] += KERNEL_VALUE(this,category,i,j);

            total_integral += KERNEL_TOTAL_VALUE(this,i,j);
        }
    }

    printf("\nKernel Weight Sums:\n");
    printf("Primary   First     Second    Multiple  Brem+Annih  Total\n");
    for (int i=0;i<N_KERNEL_CATEGORIES;i++) {
        printf("%6.2e  ", integral[i]);
    }
    printf("%6.2e\n", total_integral);
}

// write/read HDF5
int KERNEL::writeToFile(std::string fname, bool verbose) {
    // open hdf5 file
    H5::H5File h5file = H5::H5File(fname, H5F_ACC_TRUNC);
    H5::Group rootgroup = h5file.openGroup("/");
    if (!_writeToHDF5(rootgroup)) {
        if (verbose){ std::cout << "Failed to write KERNEL to \""<<fname<<"\""<<std::endl; }
        return false;
    }
    return true;
}
int KERNEL::readFromFile(KERNEL& kern, std::string fname, bool verbose) {
    H5::Exception::dontPrint();
    try {
        kern = KERNEL{};
        H5::H5File h5file = H5::H5File(fname, H5F_ACC_RDONLY);
        H5::Group rootgroup = h5file.openGroup("/");
        if (!KERNEL::_readFromHDF5(kern, rootgroup)) {
            if (verbose){ std::cout << "Failed to read KERNEL from \""<<fname<<"\""<<std::endl; }
            return false;
        }
    }
    catch (H5::FileIException &file_exists_error) {
        if (verbose){ std::cout << "Failed to read KERNEL from \""<<fname<<"\""<<std::endl; }
        return false;
    }
    return true;
}
int KERNEL::_readFromHDF5(KERNEL& kern, H5::Group& h5group) {
    { // read kernel_file string
        auto att = h5group.openAttribute("kernel_file");
        H5std_string buf("");
        att.read(att.getDataType(), buf);
        kern.kernel_file = buf;
    }
    { // read angles, radii
        auto att = h5group.openAttribute("nradii");
        att.read(H5::PredType::NATIVE_INT, &kern.nradii);
        att = h5group.openAttribute("ntheta");
        att.read(H5::PredType::NATIVE_INT, &kern.ntheta);

        att = h5group.openAttribute("radial_boundaries_cm");
        att.read(H5::PredType::NATIVE_FLOAT, &kern.radial_boundary);

        att = h5group.openAttribute("angular_boundaries_deg");
        att.read(H5::PredType::NATIVE_FLOAT, &kern.angular_boundary);
    }
    { // read matrix
        auto dset = h5group.openDataSet("weights");
        auto file_space = dset.getSpace();
        hsize_t dims[3] {};
        file_space.getSimpleExtentDims(dims);
        if (dims[0] != N_KERNEL_CATEGORIES || dims[1] != kern.ntheta || dims[2] != kern.nradii) {
            throw std::runtime_error("dataset dimensions do not agree with kernel dimensions");
        }
        hsize_t subcount[] = {1, uint(kern.ntheta), uint(kern.nradii)};
        auto mem_space = H5::DataSpace(3, subcount);
        for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
            hsize_t offset[] = {uint(cc), 0, 0};
            kern.matrix[cc] = std::vector<float>(kern.nradii * kern.ntheta);
            file_space.selectHyperslab(H5S_SELECT_SET, subcount, offset);
            dset.read(kern.matrix[cc].data(), H5::PredType::NATIVE_FLOAT, mem_space, file_space);
        }
    }
    { // read total matrix
        auto dset = h5group.openDataSet("total_weights");
        auto file_space = dset.getSpace();
        hsize_t dims[2] {};
        file_space.getSimpleExtentDims(dims);
        if (dims[0] != kern.ntheta || dims[1] != kern.nradii) {
            throw std::runtime_error("dataset dimensions do not agree with kernel dimensions");
        }
        for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
            kern.total_matrix = std::vector<float>(kern.nradii * kern.ntheta);
            dset.read(kern.total_matrix.data(), H5::PredType::NATIVE_FLOAT, file_space, file_space);
        }
    }

    return true;
}
int KERNEL::_writeToHDF5(H5::Group& h5group) const {
    H5::DataSpace scalarspace {};

    { // write kernel_file as string attr
        std::string kernel_file = std::string(this->kernel_file);
        H5::StrType str_t(H5::PredType::C_S1, kernel_file.size()+1);
        auto att = h5group.createAttribute("kernel_file", str_t, scalarspace);
        att.write(str_t, kernel_file);
    }
    { // write angles, radii
        auto att = h5group.createAttribute("nradii", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_INT, &nradii);
        att = h5group.createAttribute("ntheta", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_INT, &ntheta);

        hsize_t dims[] = { uint(nradii) };
        H5::DataSpace vectspace(1, dims);
        att = h5group.createAttribute("radial_boundaries_cm", H5::PredType::IEEE_F32LE, vectspace);
        att.write(H5::PredType::NATIVE_FLOAT, &radial_boundary);

        dims[0] = uint(ntheta);
        vectspace = H5::DataSpace(1, dims);
        att = h5group.createAttribute("angular_boundaries_deg", H5::PredType::IEEE_F32LE, vectspace);
        att.write(H5::PredType::NATIVE_FLOAT, &angular_boundary);
    }
    { // write matrix
        hsize_t dims[] = {N_KERNEL_CATEGORIES, uint(ntheta), uint(nradii)};
        H5::DataSpace file_space(3, dims);
        auto dset = h5group.createDataSet("weights", H5::PredType::IEEE_F32LE, file_space);
        hsize_t subcount[] = {1, uint(ntheta), uint(nradii)};
        H5::DataSpace mem_space(3, subcount);
        for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
            hsize_t offset[] = {uint(cc), 0, 0};
            file_space.selectHyperslab(H5S_SELECT_SET, subcount, offset);
            dset.write(matrix[cc].data(), H5::PredType::NATIVE_FLOAT, mem_space, file_space);
        }
    }
    { // write total matrix
        hsize_t dims[] = {uint(ntheta), uint(nradii)};
        H5::DataSpace dspace(2, dims);
        auto dset = h5group.createDataSet("total_weights", H5::PredType::IEEE_F32LE, dspace);
        dset.write(total_matrix.data(), H5::PredType::NATIVE_FLOAT, dspace, dspace);
    }
    return true;
}
void KERNEL::normalize(bool verbose) {
    double normkern=0.f;
    //for each kernel element, read values for all categories
    //at each category increment total kernel value for voxel i,j
    for (int tt=0; tt<ntheta; tt++) {
        for (int rr=0; rr<nradii; rr++) {
            normkern += KERNEL_TOTAL_VALUE(this, rr, tt);
        }
    }
    //
    // normalize kernel - ensure sum == 1 for total_matrix
    double normcheck = 0.0;
    for (int t=0; t<ntheta; ++t) {
        for (int r=0; r<nradii; ++r) {
            for (int cc=0; cc<N_KERNEL_CATEGORIES; cc++) {
                KERNEL_VALUE(this, cc, r, t) /= (float)normkern;
            }
            KERNEL_TOTAL_VALUE(this,r,t) /= (float)normkern;
            normcheck += KERNEL_TOTAL_VALUE(this,r,t);
        }
    }
    if (verbose) {
        printf("Pre-normalization kernel weight sum:  %g\n", normkern);
        printf("Post-normalization kernel weight sum: %g\nn", normcheck);
    }
    if (fabs(normcheck-1.f) >= 1e-6f) {
        throw std::runtime_error("normalization of kernel failed");
    }
}

