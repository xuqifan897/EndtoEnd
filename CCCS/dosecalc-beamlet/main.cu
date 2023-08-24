/* This file contains Cuda runtime calls (cuda*****()) and must be compiled directly using nvcc.
 *  Extension .cu ensures this is handled by cmake automatically
*/

// header file includes
#include <cstring> // memset()
#include <cstdio> // fopen, FILE...
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <string>
#include <exception>
#include <dirent.h>
#include <thread>

#include "H5Cpp.h" // Namespace H5::

// cuda multithreading helpers
#include "Utilities/multithreading.h"

// basic utility classes
#include "Utilities/logging.h"
#include "Utilities/timing.h"

#include "RTClasses/rtstruct.h"
#include "CudaUtilities/manage_gpu.cuh"
#include "dosecalc_defs.h"
#include "DoseCalcIO/dosecalcio.h"
#include "sparsify_manager.h"
#include "version.h"
#include "profile.h"

// kernel code is organized using standard c++ code/header organization. *_host.cu contains host program
// callable functions and data structures that must be accessible directly by host code (this module).
// *_host.cuh contains function and data declarations to be included by host code
// _device.cu/.cuh contain kernels and kernel forward declarations as well as data structures that should be
// accessible only from device code and calling functions in *_host.cu
#include "nvbbRayConvolve_host.cu"
#include "server/brain_defs.h"  // dosecalc data container structs

// GLOBAL VARS
bool extra_verbose = false; // print findBEV verbose output
bool extra_debug = false;    // print between-kernel-execution data structures (dose, terma, density)
bool TERMA_ONLY = false; // XXX temporary

static float            sparsity_thresh = 1e-6;
SparsifyManager         sparsifymanager;
ROIMaskList             roi_list;

CONSTANTS*   constants;
BEAM* device_beam_arr[MAXIMUM_DEVICE_COUNT];
MONO_KERNELS mono_kernels;

void* tworker_radConvolveTexture(void* args) {
    // cast args to struct
    DEVICE_THREAD_DATA* tdata = (DEVICE_THREAD_DATA *)(args);

    // set device number
    if (tdata->verbose && extra_verbose) {
        printf("New thread spawned, setting device to: %d\n", tdata->gpuid);
    }
    checkCudaErrors(cudaSetDevice(tdata->gpuid));

    /////////////////////// RUN CONVOLUTIONS
    // each thread executes its device in parallel
    // each thread launches streams (gpu workers) for each beam assigned to its device to maximize
    // device usage by loading workers from pool of work
    // int threadid = 0; // DEBUG
    radconvolveTexture(
            &mono_kernels,
            constants,
            tdata->device_beam_arr,
            tdata->device_nbeams,
            tdata->nrays,
            tdata->deviceid,
            tdata->gpuid,
            tdata->n_unpackstreams,
            tdata->verbose,
            tdata->timing,
            tdata->debugwrite
            );
    // wait for all processes on this GPU to complete
    checkCudaErrors( cudaDeviceSynchronize() );
    return nullptr;
}

int main(int argc, char *argv[])
{
    Logger logger; // standard stdout printing
    AutoTimer timer_total;
    AutoTimer timer_task;
    NVTX_PUSH_RANGE("0-Full Execution", 0);
    NVTX_PUSH_RANGE("1-Preparation", 1);

	// flags and variable definitions
    bool verbose = false;
    bool timing  = false;
    bool debugwrite = false;

    int num_threads_total = std::thread::hardware_concurrency();

	// output usage if prompted by -help (-h) command line flag
    if (checkCmdLineFlag( argc, (const char**)argv, "help" )) {
        printf("dosecalc-beamlet (v%s):\n"
               "  a beamlet-based, general purpose, multi-GPU dose calculation engine designed\n"
               "  to provide efficient evaluation of per-beamlet dose coefficients from customizable fluence map\n"
               "  specifications.\n\n", VERSION_STRING);
        printf("Usage:  dosecalc-beamlet [options...]\n");
        printf(" Optional Args:\n");
        printf("  --out=[outfile]                  name of result file\n");
        printf("  --ndevices=[<int>]               maximum number of GPUs to employ\n");
        printf("  --device=[<int>]                 device ID for first GPU\n");
        /* printf("  --terma                          compute beamlet terma only\n"); */
        /* printf("  --sparsify-strategy              sparsify strategy; choose from:  ['inline', 'threaded']\n"); */
        /* printf("  --write-strategy                 output write strategy; choose from:  ['central', 'perbeamlet']\n"); */
        printf("  --sparsify-threshold=[<float>]   sparsify threshold [default: %0.2e]\n", sparsity_thresh);
        /* printf("  --sparsify-streams=[<int>]       # streams used to write sparse beamlet dose to file\n"); */
        printf("  --srworkers=[<int>]              # threaded workers to use during beamlet sparsify\n");
        /* printf("  --w2fworkers=[<int>]             # threaded workers to use during beamlet write-to-file\n"); */
        printf("  --verbose[-extra]                output calculation progress information\n");
        printf("  --timing                         output calculation timing information\n");
        /* printf("  --debug[-extra]                  save debug volumes\n"); */
        printf("  --help                           show this help message\n");
        exit(EXIT_FAILURE);
    }

	// check and set flags
    /* if (checkCmdLineFlag( argc, (const char**)argv, "terma" )) { */
    /*     TERMA_ONLY = true; */
    /* } */
    if (checkCmdLineFlag( argc, (const char**)argv, "verbose" )) {
        verbose = true;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "verbose-extra" )) {
        verbose = extra_verbose = true;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "timing" )) {
        timing = true;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "debug" )) {
        debugwrite = true;
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "debug-extra" )) {
        debugwrite = extra_debug = true;
    }

///////////////////////// LOAD PRE COMPUTED DATA (from omni-precomp)
    // holds beam delivery / convolution / data volume information
    constants = new CONSTANTS{};

    // holds 3D data volumes for convolution calculation (terma, density, kernels)
    SHM_DATA* datavols = new SHM_DATA{};

    logger.print_head("LOADING DATA");

    // this data is generated by the mgcs-omni-precomp script and stored locally
    // will need to be updated to calculate on demand per beam geometry
    if ( load_omni_header( constants, verbose ) < 0) {
        printf("\n Failed to load header data.\n");
        delete constants;
        return -1;
    }
    if ( load_data( constants, datavols, verbose ) < 0) {
        printf("\n Failed to mmap pre-computed data.\n");
        delete constants;
        return -1;
    }

//////////////////////////// PARSE INPUTS
    // allow cli selection of output filename in the results folder
    std::ostringstream full_h5_fname;
    {
        std::string outfile("Dose_Coefficients");
        char* c_outfile;
        if (getCmdLineArgumentString(argc, (const char**)argv, "out", &c_outfile)) {
            outfile = std::string(c_outfile);
        }
        if (dcio::is_unqualified(outfile)) {
            if (!dcio::dir_exists(RESULTS_DIR)) {
                dcio::create_directory(RESULTS_DIR);
            }
            full_h5_fname << RESULTS_DIR << "/" << outfile << ".h5";
        } else { full_h5_fname << outfile << ".h5"; }
    }

    int n_unpackstreams = 4;
    if (checkCmdLineFlag( argc, (const char**)argv, "sparsify-streams" )) {
        n_unpackstreams = getCmdLineArgumentInt( argc, (const char**)argv, "sparsify-streams");
    }
    n_unpackstreams = max(1, n_unpackstreams);

    int nrays = constants->nphi * (constants->ntheta / 2);  // convolution directions per beamlet

    // number of GPUs to be utilized by this server
    int ndev_requested = DEFAULT_DEVICE_COUNT;
    if (checkCmdLineFlag( argc, (const char**)argv, "ndevices" )) {
        ndev_requested = getCmdLineArgumentInt( argc, (const char**)argv, "ndevices");
    }

	// sets the device ID for the first GPU
	// when using multiple GPUs, the subsequent device IDs will be incremented from this
    int DEVICE = 0;
    if (checkCmdLineFlag( argc, (const char**)argv, "device" )) {
        DEVICE = getCmdLineArgumentInt( argc, (const char**)argv, "device");
    }

    SPARSIFY_STRATEGY sparsify_strategy;
    std::string sparsify_strat;
    WRITE_STRATEGY write_strategy;
    std::string write_strat;
    {
        char* strategy;
        if (getCmdLineArgumentString(argc, (const char**)argv, "sparsify-strategy", &strategy)) {
            sparsify_strat = std::string(strategy);
            sparsify_strat = dcio::tolower(sparsify_strat);
        } else { sparsify_strat = std::string(); }
        if (sparsify_strat.compare("inline") == 0) {
            sparsify_strategy = SPARSIFY_STRATEGY::INLINE;
            sparsify_strat = "Inline";
        } else { // if (sparsify_strat.find("threaded") != std::string::npos) {
            sparsify_strategy = SPARSIFY_STRATEGY::THREADED;
            sparsify_strat = "Threaded";
        }

        if (getCmdLineArgumentString(argc, (const char**)argv, "write-strategy", &strategy)) {
            std::cout << std::string(strategy) << std::endl;
            write_strat = std::string(strategy);
            write_strat = dcio::tolower(write_strat);
        } else { write_strat = std::string(); }
        if (write_strat.compare("perbeamlet") == 0) {
            write_strategy = WRITE_STRATEGY::PER_BEAMLET;
            write_strat = "Per-Beamlet";
        /* } else if (write_strat.compare("perbeam") == 0) { */
        /*     write_strategy = WRITE_STRATEGY::PER_BEAM; */
        /*     write_strat = "Per-Beam"; */
        } else {
            write_strategy = WRITE_STRATEGY::CENTRAL;
            write_strat = "Centralized";
        }
    }

    const char* keys_sparsity_thresh[] = {"sparsify-threshold", "thresh", "sparsity-threshold", NULL};
    for (int ii=0; keys_sparsity_thresh[ii]!=NULL; ii++) {
        const char* key = keys_sparsity_thresh[ii];
        if (checkCmdLineFlag( argc, (const char**)argv, key )) {
            sparsity_thresh = getCmdLineArgumentFloat( argc, (const char**)argv, key);
            break;
        }
    }

    int num_srworkers  = 4;
    int num_w2fworkers = 1;
    if (checkCmdLineFlag( argc, (const char**)argv, "srworkers" )) {
        num_srworkers = getCmdLineArgumentInt( argc, (const char**)argv, "srworkers");
    }
    if (checkCmdLineFlag( argc, (const char**)argv, "w2fworkers" )) {
        num_w2fworkers = getCmdLineArgumentInt( argc, (const char**)argv, "w2fworkers");
    }

	// automatically load the omni-beam-list generated by the precompution script
	// function defined in "../mgcs-brain/io_functions.h"
    BEAM* beam_arr = new BEAM[constants->beam_count];

    if ( load_omni_beam_list( beam_arr, constants->beam_count, int(verbose)+2*int(extra_verbose) ) < 0 ) {
        printf("\n Failed to load beam list.\n");
        delete [] beam_arr;
        return -1;
    }
    // limit ndevices to nbeams
    if ((int)constants->beam_count < ndev_requested) {
        ndev_requested = constants->beam_count;
    }

	// load spectrum files
    mono_kernels.spectrum_file = std::string(constants->beam_spec);

    //monoenergetic kernel file
    //spectrum data
    if ( (1 != read_spectrum_file(&mono_kernels,verbose)) ) {
        printf("Failed reading spectrum file!\n");
        exit(1);
    }
    //read all monoenergetic kernels that were specified in spectrum file
    for (int i=0; i<mono_kernels.nkernels; i++) {
        if ( (1 != read_kernel(&(mono_kernels.kernel[i]))) ) {
            printf("Failed reading kernel!\n");
            exit(1);
        }
    }

    // LOAD ROIs from file - into global roi_list
    std::ostringstream roi_list_inpath;
    roi_list_inpath << Paths::Instance()->temp_dir() << "/" << "roi_list.h5";
    if (constants->reduce) {
        std::cout << "LOADING OPTIMIZATION STRUCTURES:" << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        if (ROIMaskList::readFromFile(roi_list, roi_list_inpath.str())) {
            std::cout << "Discovered " << roi_list._coll.size() << " ROIs in \""<<roi_list_inpath.str()<<"\":" << std::endl;
            uint idx = 0;
            for (const auto& v : roi_list._coll) {
                ++idx;
                std::cout << "  " << idx << ": "<< v->name << std::endl;
            }
        }
        std::cout << "--------------------------------------" << std::endl;
        std::cout << std::endl;
    }


    if (timing) {
        timer_task.restart_print_time_elapsed("Read_Specfile & Data Prep");
    }

    logger.print_tail();

////////////////////////////// INITIALIZE CUDA DEV ////////////////////////////////////
    logger.print_head("DEVICE SUMMARY");
    int ndevices = 0;
    int ndev = 0;
    int gpuid_arr[MAXIMUM_DEVICE_COUNT] = {0};
    /* init_devices_uva(ndev_uva, gpuid_arr, MAXIMUM_DEVICE_COUNT, ndev_requested, DEVICE, verbose); */
    init_devices(ndev, gpuid_arr, MAXIMUM_DEVICE_COUNT, ndev_requested, DEVICE, verbose);

    // ndev should always contain the ACTUAL number of devices being used for execution
    if (ndev == 1) {
        printf("Using 1 Device\n");
        ndevices = 1;
    }
    else if (ndev_requested < ndev) {
        printf("%d Devices, but only using %d as requested.\n", ndev, ndev_requested);
        ndevices = ndev_requested;
    } else {
        printf("Using %d Devices.\n", ndev);
        ndevices = ndev;
    }

    if (timing) {
        std::cout << std::endl;
        timer_task.restart_print_time_elapsed("Device initialization");
    }
    logger.print_tail();


////////////////////////////// ASSIGN BEAMS TO DEVICES ////////////////////////////////////
	// determine how many convolutions each GPU will perform
    // TODO: convert all dynamic mem allocations to Vector types with automatic destruction
    int device_nbeams[MAXIMUM_DEVICE_COUNT] = {0}; // total nbeams for all streams on this GPU
    int _assigned_beams = 0; // total nbeams assigned
    for (int i=0; i<ndevices; i++) {
        int _nbeams = constants->beam_count / ndevices;
        device_nbeams[i] = _nbeams;  //intrinsic floor() for positives
        _assigned_beams += _nbeams;
    }
    //assign leftovers, distributed among all devices (to stream 0)
    int _beams_left = constants->beam_count - _assigned_beams;
    for (int k=0; k<_beams_left; k++) {
        ++device_nbeams[k];
    }
    // store device specific beams into device specific arrays
    // note: array of pointers is used (rows are static, cols are dynamic)
    int beam_idx = 0;
    for (int i=0; i<ndevices; i++) {
        device_beam_arr[i] = new BEAM[ device_nbeams[i] ];
        for (int j=0; j<device_nbeams[i]; j++) {
            device_beam_arr[i][j] = beam_arr[beam_idx++];
        }
    }

    if (verbose){
        logger.print_head("BEAM ASSIGNMENTS");
        for (int i=0; i<ndevices; i++ ) {
            int gpuid = gpuid_arr[i];
            printf("GPU %d: has %d beams assigned\n", gpuid, device_nbeams[i]);
            for (int j=0; j<device_nbeams[i]; j++) {
                std::cout << "   beam " << (j+1) << ": " << device_beam_arr[i][j] << std::endl;
            }
            std::cout << std::endl;
        }
        logger.print_tail();
    }

////////////////////////////// INITIALIZE DEVICE MEMORY AND TEXTURES ////////////////////////////////////
    logger.print_head("MEMORY INITIALIZATION");

    // Store copies of constant problem data and texture/surface references on each GPU
    for (int i=0; i<ndevices; i++) {
        int gpuid = gpuid_arr[i];
        std::cout << "Initializing memory on GPU: " << gpuid << std::endl;
        checkCudaErrors(cudaSetDevice(gpuid));
        initCudaConstandTex(
                datavols,
                &mono_kernels,
                constants,
                nrays,
                i,
                gpuid,
                verbose,
                timing,
                debugwrite
                );
        std::cout << std::endl;
    }

    if (timing) {
        std::cout << std::endl;
        timer_task.restart_print_time_elapsed("GPU initialization");
    }
    logger.print_tail();

////////////////////////////// DOSE CONVOLUTION ////////////////////////////////////
    logger.print_head("DOSE CONVOLUTION");


    // populate patient metadata for storage later
    HEADER_PATIENT patient_header {};
    patient_header.N_beams = constants->beam_count;
    patient_header.dicom_start_cm = constants->start;
    patient_header.full_dicom_size = constants->size;
    patient_header.voxel_size_cm = constants->voxel;
    patient_header.rev_latspacing_cm = constants->rev_latspacing;
    patient_header.rev_longspacing_cm = constants->rev_longspacing;
    patient_header.bbox_start = constants->calc_bbox_start;
    patient_header.bbox_size = constants->calc_bbox_size;
    patient_header.penumbra_cm = constants->penumbra;
    patient_header.sparsity_thresh = sparsity_thresh;
    patient_header.nphi = constants->nphi;
    patient_header.ntheta = constants->ntheta;
    patient_header.nradii = constants->nradii;
    patient_header.kernel_extent_cm = constants->kernel_extent;
    patient_header.beam_spectrum = std::string(mono_kernels.spectrum_file);
    patient_header.target_structure = std::string{constants->target_structure};
    // Add only if reduction is to be used after sparsify()
    if (!roi_list._coll.empty()) {
        patient_header.roi_order = roi_list.getROINames();
        patient_header.row_block_capacities = roi_list.getROICapacities();
    }

    if (sparsify_strategy != SPARSIFY_STRATEGY::THREADED) {
        // TODO: Add InlineSparsifyWorker that processes on ->push() synchronously
        //       nvbbRayConvolve_host shouldn't need to know anything other than how to push data to a manager
        std::cerr << set_color(COLOR::RED) << "Sparsify Strategy: " << sparsify_strat << " not implemented" << set_color() << std::endl;
        throw std::exception();
    }

    int num_threads_available = num_threads_total - ndevices;
    /* int total_workers = num_srworkers + num_w2fworkers; */

    /* while (total_workers > num_threads_available) { */
    /*     --num_w2fworkers; */
    /*     --total_workers; */
    /*     if (total_workers > num_threads_available) { */
    /*         --num_srworkers; */
    /*         --total_workers; */
    /*     } */
    /* } */
    /* num_w2fworkers = max(1, num_w2fworkers); */
    num_srworkers = max(1, min(num_threads_available-num_w2fworkers, num_srworkers));
    n_unpackstreams = num_srworkers;

    sparsifymanager.init(full_h5_fname.str(), patient_header, num_srworkers, num_w2fworkers, sparsity_thresh, write_strategy);
    sparsifymanager.activate();
    // enabling compression slows down write-to-file drastically and doesnt save on memory
    /* sparsifyworker->set_compress_lvl(6); */
    /* sparsifyworker->set_chunksize(25); */

////////////////////////////// Print Details ////////////////////////////////////////////
    if (verbose) {
        logger.print_head("EXECUTION SUMMARY");
        printf("%3d NVB directions per beam\n", nrays);
        printf("%3d devices for server\n", ndevices);
        printf("%3d beams per device (max)\n", device_nbeams[0]);
        printf("%3d CPU threads in use (of %d)\n\n", num_srworkers + num_w2fworkers + ndevices, num_threads_total);

        std::cout << "Sparsify Strategy:     \"" << sparsify_strat << "\"" << std::endl;
        std::cout << "Data Write Strategy:   \""<<write_strat<<"\"" << std::endl;
        std::cout << "Sparsity Threshold:    " << std::scientific << std::setprecision(2) << sparsity_thresh << std::endl;
        std::cout << "Sparsify Workers:      " << sparsifymanager.get_num_srworkers() << std::endl;
        std::cout << "Write-to-file workers: " << sparsifymanager.get_num_w2fworkers() << std::endl;
        std::cout << "Sparsify Streams:      " << n_unpackstreams << std::endl;
        logger.print_tail();
    }
/////////////////////////////////////////////////////////////////////////////////////////

    // using CUTThread library for portability between windows and linux
    CUTThread* dev_threads = new CUTThread[ndevices];
    DEVICE_THREAD_DATA** tdata = new DEVICE_THREAD_DATA*[ndevices];
    int _threadid = 0;
    for (int i=0; i<ndevices; i++) {
        int gpuid = gpuid_arr[i];

        // initialize thread data
        tdata[_threadid] = new DEVICE_THREAD_DATA();  // zero initialize dynamically allocated struct

        tdata[_threadid]->device_beam_arr = device_beam_arr[i];
        tdata[_threadid]->device_nbeams   = device_nbeams[i];
        tdata[_threadid]->nrays           = nrays;
        tdata[_threadid]->deviceid        = i;
        tdata[_threadid]->gpuid           = gpuid;
        tdata[_threadid]->n_unpackstreams = n_unpackstreams;
        tdata[_threadid]->verbose         = verbose;
        tdata[_threadid]->timing          = timing;
        tdata[_threadid]->debugwrite      = debugwrite;

        // create thread and execute worker
        dev_threads[_threadid] = cutStartThread(tworker_radConvolveTexture, (void*) tdata[_threadid]);
        _threadid++;
    }

    NVTX_POP_RANGE;
    NVTX_PUSH_RANGE("1-Dose Calculation", 1);

    // join cpu threads
    cutWaitForThreads(dev_threads, ndevices);
    if (sparsify_strategy == SPARSIFY_STRATEGY::THREADED) {
        Timer timer_waiting;
        if (timing) { timer_waiting.start(); }
        sparsifymanager.deactivate();
        if (timing) { timer_waiting.stop_print_time_elapsed("Waiting for final data to be saved"); }
    }

    NVTX_POP_RANGE;
    NVTX_PUSH_RANGE("1-Cleanup", 1);

    // cleanup thread data
    delete [] dev_threads;
    for (int i=0; i<ndevices; i++) {
        delete tdata[i];
    }
    delete [] tdata;

    if (timing) {
        std::cout << std::endl;
        timer_task.restart_print_time_elapsed("Convolution");
    }
    logger.print_tail();

    logger.print_head("Post-Processing & Cleanup");
    /*******************************************************************************************************/

    // print output target
    std::cout << "Data written to \""<<full_h5_fname.str()<<"\""<<std::endl;

    // move mask file if exists
    if (dcio::file_exists(roi_list_inpath.str())) {
        std::ostringstream roi_list_outpath;
        roi_list_outpath << dcio::get_dirname(full_h5_fname.str()) << "/" << dcio::splitext(full_h5_fname.str())[0] << ".mask";
        if (dcio::copy_file(roi_list_inpath.str(), roi_list_outpath.str()))
        { std::cout << "Masks written to: "; }
        else { std::cout << "Failed writing masks to:"; }
        std::cout << "\""<< roi_list_outpath.str() <<"\"" << std::endl;
    }

    // clean up GPU(s)
    freeCudaTexture(ndevices, gpuid_arr, debugwrite);
    for (int deviceid=0; deviceid<ndevices; deviceid++) {
        checkCudaErrors( cudaSetDevice(gpuid_arr[deviceid]) );
        checkCudaErrors( cudaDeviceReset() );
    }
    free_data(constants, datavols); // unmap files from memory
    for (int deviceid=0; deviceid<ndevices; deviceid++) {
        delete [] device_beam_arr[deviceid];
    }
    delete [] beam_arr; // can probaly be done right after assigning beams since all structs are copied;
    delete datavols;
    delete constants;

    logger.print_tail();
    if (timing) {
        timer_total.stop_print_time_elapsed("Full Program Execution");
    }

    NVTX_POP_RANGE;
    NVTX_POP_RANGE;

    return 0;
}
