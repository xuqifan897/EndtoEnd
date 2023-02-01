#ifndef __ARGPARSER_H__
#define __ARGPARSER_H__

#include <string>
#include <cstring>
#include <iostream>
#include <cstdio>

#include "dosecalc_defs.h"
#include "DoseCalcIO/dosecalcio.h"
#include "Utilities/logging.h"
#include "Utilities/configparser.h"
#include "version.h"

#include "helper_cuda.h"
#include "helper_math.h"

struct CLIArgs {
    /* Quality Settings */
    std::string  config_file;
    unsigned int ntheta              = 8;             // number of azimuthal conv. angles (kernel fan)
    unsigned int nphi                = 8;             // number of zenithal conv.angles (kernel rotation)
    float        voxelsize           = 0.25f;         // isotropic voxelsize of output [unit: cm]
    float        rev_latspacing      = voxelsize;     // lateral voxel size assumed during convolution in Rays Eye View coord system, lateral spacing of convolution rays
    float        rev_longspacing     = voxelsize;     // convolution step size along convolution ray [unit: cm]
    float        penumbra            = 1.f;           // expansion of dose subvolume computation by PENUM in each of 4 directions (+/- x,z BEV)  [unit: cm]
    float        kernel_extent       = 2.f;           // dose kernel radius truncation parameter [unit: cm]
    uint         ss_factor           = 3;             // terma anti-aliasing super-sampling factor (must be >=1; 1 to disable)
    uint3        max_rev_size        = {500,500,500}; // size of REV dynamic allocation on GPU

    /* Problem Description */
    std::string  dicom_dir;
    std::string  densvol;
    std::string  beam_file;                           // path to beam_list.txt file describing beams to calculate
    int          nbeams_requested    = -1;            // use all available beams by default
    std::string  ptv_name;
    bool         target_exact_match  = false;         // match by substring by default
    std::string  bbox_roi_name;
    std::string  beam_spec           = "spec_6mv";
    std::string  ctlut_file;                          // path to file defining CT number to mass density conversion

    /* Beamlet Params */
    std::string  structures_file;                     // path to structures.json describing optimization structures
    float2       beamlet_size        = {0.5f, 0.5f};  // bi-directional size of beamlet [unit: cm]
    uint2        fmap_dims           = {40,40};       // bi-directional size of fluence map

    /* Full Beam Params */
    std::string  fmaps_file;                          // path to file containing beam metadata and manual fluence maps

    /* Hardware Settings */
    unsigned int dev                 = 0;             // GPU id

    /* Flags */
    bool         reduce              = false;         // reduction prunes all voxels from sparse dose that aren't in specified structures
    bool         apertureready       = false;         // perform aperture-ready post-processing on automatically defined fluence maps
    bool         all_beamlets        = false;         // activate all beamlets in field rather than just those touching target
    bool         nolut               = false;         // disable use of the CT LUT (hotfix)
    bool         timing              = false;
    int          verbose             = false;
    bool         debug               = false;
};

void print_usage() {
    // printf("Usage %s:\n  -phant=[CT file base name] (built-in phantom choices: water, mediastinum, slab)\n",argv[0]);
    printf("dosecalc-preprocess (v%s):\n"
           "  Preprocessor to the full/beamlet-based dose calculation strategies that handles dicom data interpretation,\n"
           "  density volume conversion and resampling, ROI mask generation, quality setting, and beam specification\n"
           "\n", VERSION_STRING);
    printf("Usage:  dosecalc-preprocess [args...]\n");
    printf(" Required Args:\n");
    printf(" [one of]:\n");
    printf("    --dicom=<dir>                    path to directory containing .dcm slices and .dcm rtstruct files\n");
    printf("    --density=<file>                 path to density volume file written in H5 format (requires use of --fmaps)\n");
    printf(" [one of]:\n");
    printf("    --beamlist=<file>                file listing beams angles\n");
    printf("    --fmaps=<file>                   file describing beams and fluence map intensities\n");
    printf("\n");
    printf(" Optional Args:\n");
    printf("   --config=<file>                   path to optional config file with setting defs\n");
    printf("\n");
    printf("   Quality Settings:\n");
    printf("     --nphi=<int>                    number of azimuthal convolution angles (kernel rotations)\n");
    printf("     --ntheta=<int>                  number of zenithal convolution angles (collapsed cones per kernel)\n");
    printf("     --penum=<float>                 lateral expansion of dose calculation volume around beam [cm]\n");
    printf("     --voxsize=<float>               isotropic voxel size used in data resampling [cm]\n");
    printf("     --convlat=<float>               lateral spacing of convolution rays [cm] (default=voxsize)\n");
    printf("     --convstep=<float>              longitudinal convolution step size [cm] (default=voxsize)\n");
    printf("     --kernel-extent=<float>         lateral radius of beamlet context [cm]\n");
    printf("     --ssfactor=<int>                terma anti-aliasing factor (must be >=1; 1 to disable)\n");
    printf("\n");
    printf("   Data Definition:\n");
    printf("     --structures=<file>             json file listing optimization/target structures\n");
    printf("     --reduce                        reduce M matrix to A matrix using OAR structures\n");
    printf("     --spec=<file>                   name of file describing beam energy spectrum (def: \"spec_6mv\")\n");
    printf("     --target=<str>                  substring of name of PTV defined in rtstruct file\n");
    printf("     --target-exact=<str>            exact name of PTV defined in rtstruct file\n");
    printf("     --bbox-roi=<str>                reduces dose calc to the box bounding the selected ROI\n");
    printf("     --nbeams=<int>                  set the maximum number of beams\n");
    printf("     --ctlut=<file>                  file (lookup table) defining CT# to density conversion\n");
    printf("     --device=<int>                  set the GPU device\n");
    printf("\n");
    printf("   Flags:\n");
    printf("     --aperture-ready                Apply \"aperture-ready\" post-process to fluence maps\n");
    printf("     --all-beamlets                  activate all beamlets instead of only those intersecting the target\n");
    printf("     --nolut                         disable CT HU to density LUT and use linear approx instead\n");
    printf("     --verbose[-extra]               output calculation progress information \n");
    printf("     --timing                        output calculation args.timing information \n");
    // printf("     --debug                         enable debug settings and save debug volumes\n");
    printf("     --help                          show this help message\n");
    exit(EXIT_FAILURE);
}
int parse_flags(CLIArgs& args, int argc, char *argv[]) {
    // parse only for boolean flags
    if (checkCmdLineFlag( argc, (const char**)argv, "reduce"  )) { args.reduce = true; };
    if (checkCmdLineFlag( argc, (const char**)argv, "aperture-ready"  )) { args.apertureready = true; };
    if (checkCmdLineFlag( argc, (const char**)argv, "all-beamlets"  )) { args.all_beamlets = true; };
    if (checkCmdLineFlag( argc, (const char**)argv, "nolut"  )) { args.nolut = true; };

    if (checkCmdLineFlag( argc, (const char**)argv, "verbose"   )) { args.verbose = 1;  }
    if (checkCmdLineFlag( argc, (const char**)argv, "verbose-extra" )) { args.verbose = 2;  }
    if (checkCmdLineFlag( argc, (const char**)argv, "noverbose" )) { args.verbose = 0; }

    if (checkCmdLineFlag( argc, (const char**)argv, "timing"    )) { args.timing  = true;  }
    if (checkCmdLineFlag( argc, (const char**)argv, "notiming"  )) { args.timing  = false; }

    if (checkCmdLineFlag( argc, (const char**)argv, "debug"     )) { args.debug   = true;  }
    if (checkCmdLineFlag( argc, (const char**)argv, "nodebug"   )) { args.debug   = false; }
    return true;
}
int parse_quality(CLIArgs& args, int argc, char *argv[]) {
    // parse quality related args
    if (checkCmdLineFlag( argc, (const char**)argv, "nphi" )) {
        args.nphi = getCmdLineArgumentInt( argc, (const char**)argv, "nphi");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "ntheta" )) {
        args.ntheta = getCmdLineArgumentInt( argc, (const char**)argv, "ntheta");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "penum" )) {
        args.penumbra = getCmdLineArgumentFloat( argc, (const char**)argv, "penum");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "voxsize" )) {
        args.voxelsize = getCmdLineArgumentFloat( argc, (const char**)argv, "voxsize");
    } else if (checkCmdLineFlag( argc, (const char**)argv, "voxelsize" )) {
        args.voxelsize = getCmdLineArgumentFloat( argc, (const char**)argv, "voxelsize");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "convlat" )) {
        args.rev_latspacing = getCmdLineArgumentFloat( argc, (const char**)argv, "convlat");
    } else {
        args.rev_latspacing = args.voxelsize;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "convstep" )) {
        args.rev_longspacing = getCmdLineArgumentFloat( argc, (const char**)argv, "convstep");
    } else {
        args.rev_latspacing = args.voxelsize;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "kernel-extent" )) {
        args.kernel_extent = getCmdLineArgumentFloat( argc, (const char**)argv, "kernel-extent");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "ssfactor" )) {
        args.ss_factor = getCmdLineArgumentFloat( argc, (const char**)argv, "ssfactor");
    }
    return true;
}
int parse_data(CLIArgs& args, int argc, char *argv[]) {
    // parse data inputs
    char *arg;

    // MANDATORY
    if (getCmdLineArgumentString( argc, (const char**)argv, "dicom", &arg )) {
        args.dicom_dir = std::string(arg);
    }
    if (getCmdLineArgumentString( argc, (const char**)argv, "density", &arg )) {
        args.densvol = std::string(arg);
    }
    if (args.dicom_dir.empty() && args.densvol.empty()) {
        fprintf(stderr, "ERROR: No Source Data Specified. Aborting. \n");
        exit(EXIT_FAILURE);
    }
    if (!args.dicom_dir.empty() && !args.densvol.empty()) {
        args.dicom_dir.clear();
        std::cout << "both \"--dicom\" and \"--density\" were specified. Using \"--density\"." << std::endl
            << "To override, exclusively set the preferred method." << std::endl;
    }

    if (getCmdLineArgumentString( argc, (const char**)argv, "beamlist", &arg )) {
        args.beam_file = std::string(arg);
    }
    if (getCmdLineArgumentString( argc, (const char**)argv, "fmaps", &arg )) {
        args.fmaps_file = std::string(arg);
        //TODO: if using fmaps, get all quality settings from file and allow cli args to override if present
        // for now, only beam spec and fluence maps are extracted
    }
    if (args.beam_file.empty() && args.fmaps_file.empty()) {
        std::cerr << "ERROR: one of the following arguments must be specified:" << std::endl;
        std::cerr << "  --beamlist=<file>" << std::endl;
        std::cerr << "  --fmaps=<file>" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!args.beam_file.empty() && !args.fmaps_file.empty()) {
        args.beam_file.clear();
        std::cout << "both \"--beamlist\" and \"--fmaps\" were specified. Using \"--fmaps\"." << std::endl
            << "To override, exclusively set the preferred method." << std::endl;
    }

    // OPTIONAL
    if (getCmdLineArgumentString( argc, (const char**)argv, "spec", &arg )) {
        args.beam_spec = std::string(arg);
    }

    if (getCmdLineArgumentString( argc, (const char**)argv, "structures", &arg )) {
        args.structures_file = std::string(arg);
    }
    if (getCmdLineArgumentString( argc, (const char**)argv, "ctlut", &arg )) {
        args.ctlut_file = std::string(arg);
    }

    if (checkCmdLineFlag( argc, (const char**)argv, "nbeams" )) {
        args.nbeams_requested = getCmdLineArgumentInt( argc, (const char**)argv, "nbeams");
    }

    if (checkCmdLineFlag( argc, (const char**)argv, "device" )) {
        args.dev = getCmdLineArgumentInt( argc, (const char**)argv, "device");
    }

    // read user preference for target structure name
    /* If exact is used, dont consider substring matching, otherwise set default sub or let user pick */
    if (getCmdLineArgumentString( argc, (const char**)argv, "target-exact", &arg )) {
        args.ptv_name = std::string(arg);
        args.target_exact_match= true;
    } else if (getCmdLineArgumentString( argc, (const char**)argv, "target", &arg )) {
        args.ptv_name = std::string(arg);
        args.target_exact_match = false;
    }

    if (getCmdLineArgumentString( argc, (const char**)argv, "bbox-roi", &arg )) {
        args.bbox_roi_name = std::string(arg);
    }

    return true;
}
int parse_config(CLIArgs& args, int argc, char *argv[]) {
    // parse settings in config file with path given by --config
    char *arg;

    if (getCmdLineArgumentString( argc, (const char**)argv, "config", &arg )) {
        args.config_file = std::string(arg);
    }

    if (args.config_file.empty()) { return true; }
    rapidjson::Document conf;
    if (!(conf = read_config(args.config_file)).IsObject()) {
        std::cout << "Error encountered attempting to read from \""<<args.config_file<<"\"" << std::endl;
        return false;
    }

    // BEGIN TESTING VALUES
    rapidjson::Value::ConstMemberIterator itr;

    int verbose = args.verbose;
    if (verbose) { std::cout << "READING CONFIG:" << std::endl << "---------------" << std::endl; }
    get_and_assign_scalar<float>(args.voxelsize, conf, "voxsize", verbose);
    get_and_assign_scalar<float>(args.rev_latspacing, conf, "convlat", verbose);
    get_and_assign_scalar<float>(args.rev_longspacing, conf, "convstep", verbose);
    get_and_assign_scalar<uint>(args.ss_factor, conf, "ssfactor", verbose);
    get_and_assign_vector<float, float2, 2>(args.beamlet_size, conf, "beamlet-size", verbose);
    get_and_assign_vector<uint, uint2, 2>(args.fmap_dims, conf, "fmap-dims", verbose);
    get_and_assign_scalar<float>(args.penumbra, conf, "penum", verbose);
    get_and_assign_scalar<float>(args.kernel_extent, conf, "kernel-extent", verbose);
    get_and_assign_scalar<uint>(args.ntheta, conf, "ntheta", verbose);
    get_and_assign_scalar<uint>(args.nphi, conf, "nphi", verbose);
    get_and_assign_vector<uint, uint3, 3>(args.max_rev_size, conf, "max-rev-size", verbose);
    get_and_assign_scalar<int>(args.verbose, conf, "verbose", verbose);
    get_and_assign_scalar<bool>(args.timing, conf, "timing", verbose);
    get_and_assign_scalar<int>(args.nbeams_requested, conf, "nbeams", verbose);
    get_and_assign_scalar<bool>(args.reduce, conf, "reduce", verbose);
    get_and_assign_scalar<bool>(args.apertureready, conf, "apertureready", verbose);

    get_and_assign_scalar<std::string>(args.beam_spec, conf, "spec", verbose);
    if (get_and_assign_scalar<std::string>(args.ptv_name, conf, "target-exact", verbose)) {
        args.target_exact_match = true;
    } else {
        get_and_assign_scalar<std::string>(args.ptv_name, conf, "target", verbose);
    }
    get_and_assign_scalar<std::string>(args.bbox_roi_name, conf, "bbox-roi", verbose);
    get_and_assign_scalar<std::string>(args.beam_file, conf, "beamlist", verbose);
    get_and_assign_scalar<std::string>(args.structures_file, conf, "structures", verbose);
    get_and_assign_scalar<std::string>(args.ctlut_file, conf, "ctlut", verbose);
    get_and_assign_scalar<std::string>(args.fmaps_file, conf, "fmaps", verbose);
    if (verbose) { std::cout << std::endl; }

    return true;
}

int parse_args(CLIArgs& args, int argc, char *argv[]) {
    if (argc < 2 || checkCmdLineFlag( argc, (const char**)argv, "help" )) { print_usage(); }

    if (!parse_flags(args, argc, argv)) { return false; }
    if (!parse_config(args, argc, argv)) { return false; }
    if (!parse_flags(args, argc, argv)) { return false; }
    if (!parse_quality(args, argc, argv)) { return false; }
    if (!parse_data(args, argc, argv)) { return false; }
    return true;
}

#endif //__ARGPARSER_H__
