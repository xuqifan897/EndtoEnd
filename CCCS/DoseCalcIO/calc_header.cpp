#include "calc_header.h"

// #include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "dosecalc_defs.h"
#include "server/brain_defs.h"

int write_omni_header(
        uint3 count, float3 inc, float3 start,
        uint3 calc_bbox_start, uint3 calc_bbox_size,
        int nphi, int ntheta, int nradii,
        int beam_count, const std::string& beam_spec,
        const std::string& target_structure,
        float rev_latspacing,
        float rev_longspacing,
        float kernel_extent,
        uint  ss_factor,
        uint3 max_rev_size,
        float penumbra,
        bool reduce
) {
    char headerfile[1024];
    FILE *hout;

    sprintf(headerfile,"%s/omni-header.txt",Paths::Instance()->temp_dir().c_str());
    printf("\nheader file name is %s",headerfile);

    if ( (hout = fopen(headerfile,"w")) == NULL)
    {
        printf("Cannot open text file\n");
        return(0);
    }

    fprintf(hout,
            "# ====================================================================\n"
            "# ==                         OMNI-HEADER                            ==\n"
            "# ====================================================================\n"
            "# | Nx Ny Nz           - Dicom Volume Dimensions                     |\n"
            "# | dx dy dz           - Voxel Size (cm)                             |\n"
            "# |  x  y  z           - Dicom Volume Start Coords (cm)              |\n"
            "# |  i  j  k           - Dose Bounding Box Start Indices             |\n"
            "# | Bx By Bz           - Dose Bounding Box Dimensions                |\n"
            "# | Rx Ry Rz           - REV Convolution Array Dimensions            |\n"
            "# | convlat            - Convolution ray lateral spacing (cm)        |\n"
            "# | convstep           - Convolution step spacing (cm)               |\n"
            "# | kernel_extent      - Dose kernel radius truncate distance (cm)   |\n"
            "# | ss_factor          - Terma anti-aliasing (super-sampling) factor |\n"
            "# | nphi ntheta nradii - CCK Kernel Dimensions                       |\n"
            "# | penumbra           - beamlet transverse dose spread (cm)         |\n"
            "# | beam_count         - Number of beams to pick from beam_list.txt  |\n"
            "# | beam_spectrum      - beam energy spectrum file to use            |\n"
            "# | target_structure   - Name of selected target contour             |\n"
            "# | reduce coeff. mat  - M-matrix to A-matrix reduction requested?   |\n"
            "# ===================================================================\n\n"
            );

    fprintf(hout,"%u %u %u\n", count.x, count.y, count.z );
    fprintf(hout,"%f %f %f\n", inc.x, inc.y, inc.z );
    fprintf(hout,"%f %f %f\n", start.x, start.y, start.z );
    fprintf(hout,"%u %u %u\n", calc_bbox_start.x, calc_bbox_start.y, calc_bbox_start.z);
    fprintf(hout,"%u %u %u\n", calc_bbox_size.x, calc_bbox_size.y, calc_bbox_size.z);
    fprintf(hout,"%u %u %u\n", max_rev_size.x, max_rev_size.y, max_rev_size.z);
    fprintf(hout,"%f\n",       rev_latspacing );
    fprintf(hout,"%f\n",       rev_longspacing );
    fprintf(hout,"%f\n",       kernel_extent );
    fprintf(hout,"%u\n",       ss_factor );
    fprintf(hout,"%d %d %d\n", nphi, ntheta, nradii );
    fprintf(hout,"%f\n",       penumbra );
    fprintf(hout,"%d\n",       beam_count );
    fprintf(hout,"%s\n",       beam_spec.c_str());
    fprintf(hout,"%s\n",       target_structure.c_str());
    fprintf(hout,"%d\n",       (int)reduce);
    fclose(hout);

    return 1;
}
int load_omni_header( CONSTANTS *host, bool verbose )
{
    char headerfile[1024];
    sprintf(headerfile,"%s/omni-header.txt",Paths::Instance()->temp_dir().c_str());
    FILE *header_in;
    if ( (header_in = fopen(headerfile,"r")) == NULL ) {
        printf("Cannot open data header file!\n");
        return -1;
    }
    // skip header comment lines
    int linenum = 0;
    char dump[150];
    while (true) {
        linenum++;
        fgets(dump, sizeof(dump), header_in);
        if (*dump != '#' && *dump != ' ') { break; }
    }

    fscanf(header_in,"%u %u %u", &host->size.x, &host->size.y, &host->size.z );
    fscanf(header_in,"%f %f %f", &host->voxel.x, &host->voxel.y, &host->voxel.z );
    fscanf(header_in,"%f %f %f", &host->start.x, &host->start.y, &host->start.z );
    fscanf(header_in,"%u %u %u", &host->calc_bbox_start.x, &host->calc_bbox_start.y, &host->calc_bbox_start.z);
    fscanf(header_in,"%u %u %u", &host->calc_bbox_size.x, &host->calc_bbox_size.y, &host->calc_bbox_size.z);
    fscanf(header_in,"%u %u %u", &host->max_rev_size.x, &host->max_rev_size.y, &host->max_rev_size.z);
    fscanf(header_in,"%f",       &host->rev_latspacing);
    fscanf(header_in,"%f",       &host->rev_longspacing);
    fscanf(header_in,"%f",       &host->kernel_extent);
    fscanf(header_in,"%u",       &host->ss_factor);
    fscanf(header_in,"%d %d %d", &host->nphi, &host->ntheta, &host->nradii );
    fscanf(header_in,"%f",       &host->penumbra);
    fscanf(header_in,"%d",       &host->beam_count );
    fgets(dump, sizeof(dump), header_in);
    fgets(host->beam_spec,99, header_in );
    host->beam_spec[strlen(host->beam_spec)-1] = '\0';
    fgets(host->target_structure,99, header_in );
    host->target_structure[strlen(host->target_structure)-1] = '\0';
    int tempreduce;
    fscanf(header_in, "%d", &tempreduce);
    host->reduce = (bool)tempreduce;

    if (ferror(header_in) || feof(header_in)) {
        printf("Error Reading Omni-header, Aborting\n");
        exit(1);
    }
    fclose(header_in);

    // enforce valid bbox spec
    if (host->calc_bbox_start.x >= host->size.x) { host->calc_bbox_start.x = 0; }
    if (host->calc_bbox_start.y >= host->size.y) { host->calc_bbox_start.y = 0; }
    if (host->calc_bbox_start.z >= host->size.z) { host->calc_bbox_start.z = 0; }
    if (host->calc_bbox_size.x + host->calc_bbox_start.x > host->size.x) { host->calc_bbox_size.x = host->size.x - host->calc_bbox_start.x; }
    if (host->calc_bbox_size.y + host->calc_bbox_start.y > host->size.y) { host->calc_bbox_size.y = host->size.y - host->calc_bbox_start.y; }
    if (host->calc_bbox_size.z + host->calc_bbox_start.z > host->size.z) { host->calc_bbox_size.z = host->size.z - host->calc_bbox_start.z; }

    // precompute for use in CUDA kernels
    host->calc_bbox_precomp_2d_size = host->calc_bbox_size.x * host->calc_bbox_size.y;

    if (verbose) {
        printf("OMNI-HEADER DATA:\n");
        printf("                             X        Y        Z\n");
        printf("  Data Dimensions:     %7d  %7d  %7d\n",     host->size.x,  host->size.y,  host->size.z );
        printf("  Voxel Size (cm):     %7.3f  %7.3f  %7.3f\n", host->voxel.x, host->voxel.y, host->voxel.z );
        printf("  Data Start (cm):     %7.2f  %7.2f  %7.2f\n", host->start.x, host->start.y, host->start.z );
        printf("  Dose BBox Start:     %7d  %7d  %7d\n",       host->calc_bbox_start.x, host->calc_bbox_start.y, host->calc_bbox_start.z );
        printf("  Dose BBox Size:      %7d  %7d  %7d\n",       host->calc_bbox_size.x, host->calc_bbox_size.y, host->calc_bbox_size.z );
        printf("  Nphi:                   %d\n", host->nphi);
        printf("  Ntheta:                 %d\n", host->ntheta);
        printf("  Nradii:                 %d\n", host->nradii );
        printf("  Conv Lat. Spacing (cm): %0.3f\n", host->rev_latspacing);
        printf("  Conv Step Size (cm):    %0.3f\n", host->rev_longspacing);
        printf("  Kernel Cutoff (cm):     %0.3f\n", host->kernel_extent);
        printf("  Terma SS factor:        %u\n", host->ss_factor);
        printf("  Penumbra (cm):          %0.3f\n", host->penumbra);
        printf("  # Beams:                %d\n", host->beam_count );
        printf("  Beam-spec:              %s\n", host->beam_spec);
        printf("  Target Structure:       %s\n", host->target_structure);
        printf("  Reduce Dose Matrix:     %s\n\n", host->reduce ? "yes" : "no");
    }

    return 1;
}
