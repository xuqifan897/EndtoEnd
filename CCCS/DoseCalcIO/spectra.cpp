#include "spectra.h"

#include "./paths.h"
#include "./kernel.h"
#include "dosecalc_defs.h"

int read_spectrum_file(MONO_KERNELS *mono, bool verbose)
{
    int i;
    char filename[300];
    float sum_fluence;
    FILE *specfile;

    //open spectrum file
    sprintf(filename,"%s/%s.spec",Paths::Instance()->spectra_dir().c_str(),mono->spectrum_file.c_str());
    if ( (specfile = fopen(filename,"r")) == NULL){
        printf("Cannot open spectrum file\n");
        return(-1);
    }

    mono->nkernels = 0;
    sum_fluence = 0;
    if (verbose) {
        printf("SPECTRUM-DATA:\n");
        printf("  Spectrum file: %s.spec\n\n", mono->spectrum_file.c_str());
        printf("                        mu/rho\n");
        printf("  Energy  Fluence  _atten_____en__   File\n");
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
    for (i=0;i<mono->nkernels;i++)
        mono->fluence[i] /= sum_fluence;

    fclose(specfile);
    return(1);
}

