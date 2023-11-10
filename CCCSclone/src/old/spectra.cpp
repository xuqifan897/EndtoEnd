#include "spectra.h"
#include "paths.h"
#include "kernel.h"

#include <iostream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int old::read_spectrum_file(MONO_KERNELS* mono, bool verbose)
{
    std::string filename;
    float sum_fluence;
    FILE* specfile;

    filename = (fs::path(Paths::Instance()->spectra_dir()) / 
        fs::path(mono->spectrum_file)).string();
    if ( (specfile = fopen(filename.c_str(), "r")) == nullptr)
    {
        std::cout << "Cannot open spectrum file " << filename << std::endl;
        return -1;
    }

    mono->nkernels = 0;
    sum_fluence = 0.;
    if (verbose)
    {
        std::cout << "SPECTRUM-DATA:" << std::endl;
        std::cout << "  Spectrum file: " << mono->spectrum_file;
        std::cout << "                        mu/rho\n" << std::endl;
        // ex .    0.20   0.0010   0.1370   0.0297   scaf200
    }

    //read five items from each line of spectrum file
    //energy, fluence, mu, mu_en, kernel file
    char buf[300];
    while (fscanf(specfile,"%f %f %f %f %s",
                &(mono->energy[mono->nkernels]),&(mono->fluence[mono->nkernels]),
                &(mono->mu[mono->nkernels]),&(mono->mu_en[mono->nkernels]),
                buf) == 5)
    {
        mono->kernel[mono->nkernels].kernel_file = std::string(buf);
        if (verbose){
            printf("  %5.2f   %6.4f   %6.4f   %6.4f   %s\n",
                    mono->energy[mono->nkernels],mono->fluence[mono->nkernels],
                    mono->mu[mono->nkernels],mono->mu_en[mono->nkernels],mono->kernel[mono->nkernels].kernel_file.c_str());
        }

        //sum fluence
        sum_fluence += mono->fluence[mono->nkernels];

        //count kernels
        mono->nkernels++;
    }
    if (verbose) {
        printf("\n  Number of kernels: %d\n\n",mono->nkernels);
    }

    //normalise fluence to sum to unity
    for (int i=0;i<mono->nkernels;i++)
        mono->fluence[i] /= sum_fluence;

    fclose(specfile);
    return(1);
}