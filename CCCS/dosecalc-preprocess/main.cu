// Header files containing class and struct definitions
// as well as i/o functions
#include <iostream>
#include <cstring>
#include <string>

#include "dosecalc_defs.h"
#include "DoseCalcIO/dosecalcio.h"
#include "RTClasses/rtimages.h"
#include "RTClasses/rtstruct.h"
#include "CudaUtilities/manage_gpu.cuh"
#include "CudaUtilities/geometry.cuh"
#include "Utilities/logging.h"
#include "argparser.h"
#include "cudaMemFunctions.cuh"
#include "fmapProcessing.hpp"

#include "H5Cpp.h"


bool debugwrite; // write certain debug data to data/temp/* - set by cli arg --debug

void clear_temp_files() {
    // Remove leftover temp files
    std::vector<std::string> dirs = {
        Paths::Instance()->temp_dir()+"/fluence_maps/"
    };
    dcio::remove_directories(dirs, true, true);
    dcio::remove_directory(Paths::Instance()->temp_dir(), true, true);
    if (!dcio::dir_exists(Paths::Instance()->temp_dir())) {
        dcio::create_directory(Paths::Instance()->temp_dir());
    }
}

int main(int argc, char *argv[])
{
    Logger logger;
    FloatVolume       density;      // structure to hold density data
    KERNEL            poly_kernel{};  // structure to hold polyenergetic kernel data
    MONO_KERNELS      mono_kernels{}; // array of structures to hold monoenergetic kernels
    std::vector<BEAM> beams;        // vector of objects to hold field data

    clock_t start_total = clock();
    clock_t start_task = clock();


    ////////////////////* PARSE ARGS *//////////////////////
    ////////////////////////////////////////////////////////
    CLIArgs args {};                        // container for all specifyable cli args
    /* args.max_rev_size = make_uint3((unsigned int)floor(1024.0*cbrt((query_device_memory().free-3.0)/12.0))); */

    // Parse Args
    if (!parse_args(args, argc, argv)) { std::cout << "Failed to parse args" << std::endl; exit(EXIT_FAILURE); }

    // POST-PROCESS ARGS
    debugwrite = args.debug;
    mono_kernels.spectrum_file = args.beam_spec;
    setActiveGPU(args.dev);
    if (args.ntheta % 2 != 0) {
        args.ntheta += 1;
        std::cout << set_color(COLOR::YELLOW) << "ntheta must be an even number. It has been set to "<< args.ntheta <<
            " for this execution" << set_color() << std::endl;
    }
    args.ss_factor = min(8, (max(1, args.ss_factor)));
    if (!args.densvol.empty() && args.fmaps_file.empty()) {
        std::cerr << "Specification of beams, fluence maps, and beam isocenters via \"--fmaps\" file is required when using \"--densvol\"" << std::endl;
        return EXIT_FAILURE;
    }
    if ((args.nphi % 2) != 0) {
        std::cerr << set_color(COLOR::YELLOW) << "Number of requested convolution plane angles (--nphi="<<args.nphi<<") must be an even number. Setting to " <<++args.nphi<< " for this run instead" << set_color() << std::endl;
    }
    if (args.ntheta > N_KERNEL_ANGLES) {
        std::cerr << set_color(COLOR::YELLOW) << "Number of requested convolution fan angles (--ntheta="<<args.ntheta<<") exceeds the number available. Limiting to "<<N_KERNEL_ANGLES<<" instead"<< set_color() << std::endl;
        args.ntheta = N_KERNEL_ANGLES;
    }

    // Determine program flow state flags
    bool generate_fmaps;

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////

    ///////////////////* LOAD DATA *//////////////////////
    //////////////////////////////////////////////////////
    if (!args.beam_file.empty()) {
        // this variable determines the maximum number of lines (beams) to read from that file
        // if the number of beams in the file is less than this number, it should be ok
        // in this case, beam_count will be mutated to match the number in the file
        if (args.nbeams_requested > 0) { std::cout << "Requested " << args.nbeams_requested << " beams:" << std::endl; }
        else { std::cout << "Inferring number of beams from beamlist file:" << std::endl; }

        if ( !load_beam_list( beams, args.beam_file, args.nbeams_requested, args.verbose ) || beams.size() <= 0) {
            printf("Failed to load beam list. Exiting. \n");
            return -1;
        }
        std::cout << "Loaded " << beams.size() << " beams from \""<<args.beam_file<<"\"." << std::endl;
        generate_fmaps = true;
    } else if (!args.fmaps_file.empty()) {
        HEADER_PATIENT* patient_header = nullptr; // Unused
        read_fmaps(args.fmaps_file, patient_header, beams);

        for (const auto& b : beams) {
            if (args.verbose) { std::cout << "  Beam..."<<std::setw(4)<<std::setfill('.')<<b.uid<<": "<<b<< std::endl; }
        }
        std::cout << "Loaded "<< beams.size() << " beams from \""<<args.fmaps_file<<"\"." << std::endl;
        generate_fmaps = false;
    } // else is handled in argparser

    if (beams.size() <=0) {
        std::cout << "ERROR: Must provide at least one valid beam specification to proceed" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (args.ptv_name.empty()) {
        // fallback string for PTV matching - no other spec provided
        args.ptv_name = std::string("P_");
    }

    if (args.timing) { // write args.timing
        printf("\nTime elapsed during Device Query: %4.3f msec\n", ((float)clock() - start_task)*1000 / CLOCKS_PER_SEC );
        start_task = clock();
    }
    //
    // process bbox selection
    uint3 calc_bbox_start {};
    uint3 calc_bbox_size {};

    // delete temporary files and dirs
    clear_temp_files();


    if (!args.dicom_dir.empty()) {
        logger.print_head("READ DICOM");
        printf("\n Dicom Directory: %s\n",args.dicom_dir.c_str());
        RTImage rtimage = RTImage{args.dicom_dir, (bool)args.verbose};  // exception thrown upon failure
        // const std::string& patientPosition = rtimage.getSlicePatientPosition(0);
        /* if (patientPosition != "HFS" && patientPosition != "HFP") { */
        /*     std::cout << "CT FORMAT WARNING: CT data should have \"HFS\" or \"HFP\" PatientPosition, not \""<< patientPosition <<"\""<<std::endl; */
        /*     // exit(EXIT_FAILURE); */
        /* } */
        const std::array<float, 6>& imageOrientationPatient = rtimage.getSliceImageOrientationPatient(0);
        std::array<float, 6> orient_hfs = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        if (imageOrientationPatient != orient_hfs) {
            std::cout << "Warning: CT Data orientation does not match standard orientation (HFS). Please check your orientation to be sure it is what you expect: [" <<
                imageOrientationPatient[0] << "\\" <<
                imageOrientationPatient[1] << "\\" <<
                imageOrientationPatient[2] << "\\" <<
                imageOrientationPatient[3] << "\\" <<
                imageOrientationPatient[4] << "\\" <<
                imageOrientationPatient[5] << "]" <<
                std::endl;
        }

        // structure defined in "precomp_defs.h"
        FloatVolume ctdata {};
        // get dicom data size
        ctdata.size   = make_uint3(rtimage.getDataSize());
        // get voxel dimensions, convert to cm
        ctdata.voxsize     = 0.1f * rtimage.getVoxelSize();
        // get position of voxel 0,0,0 - convert to cm
        ctdata.start   = 0.1f * rtimage.getSliceImagePositionPatient( 0 );
        // HU data volume
        ctdata.set_data(rtimage.getDataArray(), ctdata.nvoxels());

        // get CT conversion lookup table
        // read in HU to mass density conversion table (as HU:dens pairs)
        auto ctLUT = CTLUT();
        if (!args.nolut) {
            ctLUT = CTLUT();
            if (!args.ctlut_file.empty()) {
                ctLUT.label = "User Specified";
                if (!load_lookup_table(ctLUT, args.ctlut_file, args.verbose)) {
                    char msg[300];
                    sprintf(msg, "Failed to read from ct lookup table: \"%s\"", args.ctlut_file.c_str());
                    throw std::runtime_error(msg);
                }
            } else {
                ctLUT.label = "Siemens (default)";
                // LUT from: http://sbcrowe.net/ct-density-tables/
                ctLUT.points.emplace_back("Air",           -969.8f, 0.f    ) ;
                ctLUT.points.emplace_back("Lung 300",      -712.9f, 0.290f ) ;
                ctLUT.points.emplace_back("Lung 450",      -536.5f, 0.450f ) ;
                ctLUT.points.emplace_back("Adipose",       -95.6f,  0.943f ) ;
                ctLUT.points.emplace_back("Breast",        -45.6f,  0.985f ) ;
                ctLUT.points.emplace_back("Water",         -5.6f,   1.000f ) ;
                ctLUT.points.emplace_back("Solid Water",   -1.9f,   1.016f ) ;
                ctLUT.points.emplace_back("Brain",         25.7f,   1.052f ) ;
                ctLUT.points.emplace_back("Liver",         65.6f,   1.089f ) ;
                ctLUT.points.emplace_back("Inner Bone",    207.5f,  1.145f ) ;
                ctLUT.points.emplace_back("B-200",         220.7f,  1.159f ) ;
                ctLUT.points.emplace_back("CB2 30%",       429.9f,  1.335f ) ;
                ctLUT.points.emplace_back("CB2 50%",       775.3f,  1.560f ) ;
                ctLUT.points.emplace_back("Cortical Bone", 1173.7f, 1.823f ) ;
            }
            ctLUT.sort();
            if (args.verbose) { std::cout << ctLUT << std::endl; }
        }

        //convert HU dicom data to float precision density values and make isotropic in size
        if ( 1 != read_dicom(density,ctdata,args.voxelsize, args.nolut?NULL:&ctLUT, args.verbose) ) {
            printf("Failed reading CT data!\n");
            exit(EXIT_FAILURE);
        }

        // store volume properties for use below
        FrameOfReference frameofref {
            density.size,
            density.start,
            density.voxsize
        };

        // Load rtstruct data
        RTStruct rtstruct {};
        rtstruct.setDicomDirectory(  args.dicom_dir.c_str() );
        if ( !rtstruct.loadDicomInfo(args.verbose) ) {
            printf(" Couldn't load rtstruct from \"%s\". exiting\n", args.dicom_dir.c_str());
            return 1;
        }

        // populate list of rois from rtstruct
        rtstruct.loadRTStructInfo(args.verbose);
        std::vector<std::string> roi_names = rtstruct.getROINames();
        logger.print_tail();

        logger.print_head("PTV SELECTION");
        int ptv_idx = getROIIndex(rtstruct, args.ptv_name, args.target_exact_match, args.verbose);
        if (ptv_idx < 0) {
            printf("No contour could be matched from search string: %s. exiting\n", args.ptv_name.c_str());
            return 1;
        } else {
            args.ptv_name = std::string{rtstruct.getROIName(ptv_idx)};
            printf("Structure found: #%d - %s\n",ptv_idx+1, args.ptv_name.c_str());
        }
        StructureSet ptv;
        if (!loadStructureSet(ptv, rtstruct, ptv_idx, args.verbose)) {
            if (!args.verbose) std::cout << "Failed to load ROI Data for: \""<< args.ptv_name <<"\"" << std::endl;
            return 1;
        }
        logger.print_tail();
        std::cout << std::endl;

        if (args.bbox_roi_name.empty()) {
            // Guess at reasonable structure to use
            // TODO: make case insensitive search
            args.bbox_roi_name = std::string("Body");
        }
        logger.print_head("CALC. BBOX SELECTION");
        bool use_default_bbox = false;
        StructureSet bbox_roi {};
        int bbox_roi_idx = getROIIndex(rtstruct, args.bbox_roi_name, true, args.verbose); //TODO Add substring matching support in CLI
        if (bbox_roi_idx < 0) {
            use_default_bbox = true;
            printf("No contour could be matched from search string: %s. Using full volume (%d, %d, %d)\n", args.bbox_roi_name.c_str(), density.size.x, density.size.y, density.size.z);
        } else {
            if (bbox_roi_idx == ptv_idx) {
                // dont need to reload data, just copy from ptv
                args.bbox_roi_name = args.ptv_name;
                printf("Reusing structure: #%d - \"%s\"\n",ptv_idx+1, args.bbox_roi_name.c_str());
                bbox_roi = ptv; // Copy via assignment operator
            } else {
                args.bbox_roi_name = std::string{rtstruct.getROIName(bbox_roi_idx)};
                printf("Structure found: #%d - \"%s\"\n",bbox_roi_idx+1, args.bbox_roi_name.c_str());
                if (!loadStructureSet(bbox_roi, rtstruct, bbox_roi_idx, args.verbose)) {
                    if (!args.verbose) std::cout << "Failed to load ROI Data for: \""<< args.bbox_roi_name <<"\"" << std::endl;
                    use_default_bbox = true;
                }
            }
        }
        if (use_default_bbox) {
            calc_bbox_start = make_uint3(0,0,0);
            calc_bbox_size = density.size;
        } else {
            // get ROI extents
            ArrayProps extents = getROIExtents(bbox_roi, frameofref);
            calc_bbox_start = extents.crop_start;
            calc_bbox_size = extents.crop_size;
        }

        // hotfix - prevent bbox from meeting x and y edges
        // TODO: This may be unnecessary since replacing mask generator with opencv function
        calc_bbox_start.x = max(calc_bbox_start.x, 1);
        calc_bbox_start.y = max(calc_bbox_start.y, 1);
        calc_bbox_start.z = max(calc_bbox_start.z, 1);
        calc_bbox_size.x = min(calc_bbox_size.x, frameofref.size.x-2);
        calc_bbox_size.y = min(calc_bbox_size.y, frameofref.size.y-2);
        calc_bbox_size.z = min(calc_bbox_size.z, frameofref.size.z-2);

        std::cout << std::endl;
        printf("-- bbox start: (%3d, %3d, %3d)\n", calc_bbox_start.x, calc_bbox_start.y, calc_bbox_start.z);
        printf("-- bbox size:  (%3d, %3d, %3d)\n", calc_bbox_size.x, calc_bbox_size.y, calc_bbox_size.z);
        logger.print_tail();
        std::cout << std::endl;

        // create conformal fluence masks
        logger.print_head("CALCULATE PTV CENTROID");
        float3 centroid {};
        Volume<uint8_t> ptv_mask;
        {
            std::cout << "Creating PTV Mask:" << std::endl;

            ptv_mask = generateContourMask(ptv, ctdata.get_frame(), density.get_frame());

            // Calculate Centroid
            centroid = getROICentroid(ptv_mask, frameofref);
            printf(" PTV Isocenter: %f  %f  %f\n",centroid.x,centroid.y,centroid.z);

            if (args.debug) {
                // Write out mask volume
                std::vector<float> fmask(ptv_mask._vect.begin(), ptv_mask._vect.end());
                write_debug_data<float>(fmask.data(), ptv_mask.size, "debug_ptv_mask", true);
            }

        }
        for (int bb=0; bb<beams.size(); bb++) {
            BEAM& beam = beams[bb];
            if (beam.isocenter_type == BEAM::ISO_T::UNSPEC) {
                // update isocenter with ptv centroid
                beam.isocenter = centroid;
                beam.isocenter_type = BEAM::ISO_T::PTV_CENTROID;
                beam.isocenter_location = BEAM::ISO_LOC_T::IN_PTV;
                // recompute source/dir with new isocenter
                beam.reconfigure();
            } else if (beam.isocenter_type == BEAM::ISO_T::MANUAL) {
                // check if isocenter is in PTV
                int3 isoidx = make_int3( (beam.isocenter - density.start)/density.voxsize );
                int idx = isoidx.x + density.size.x*(isoidx.y + density.size.y*isoidx.z);
                if ((idx>=0 && idx<density.nvoxels()) && ptv_mask[idx]!=0) {
                    beam.isocenter_location = BEAM::ISO_LOC_T::IN_PTV;
                } else {
                    beam.isocenter_location = BEAM::ISO_LOC_T::OUT_PTV;
                    if (args.verbose<2) {
                        std::cout << set_color(COLOR::YELLOW) << "Manual isocenter spec (beam #"<<bb+1<<") outside PTV" << set_color() << std::endl;
                    }
                }
                if (args.verbose>=2) {
                    if (beam.isocenter_location == BEAM::ISO_LOC_T::OUT_PTV) {
                        std::cout << set_color(COLOR::YELLOW);
                    }
                    std::cout << "  Beam..."<<std::setw(4)<<std::setfill('.')<<bb+1<<": Isocenter "<<(beam.isocenter_location==BEAM::ISO_LOC_T::IN_PTV?"INSIDE":"OUTSIDE")<<" PTV | indices=("<<isoidx.x<<", "<<isoidx.y<<", "<<isoidx.z<<")" << std::endl;
                    std::cout << set_color();
                }
            }


            if (generate_fmaps) {
                beam.fmap_size = args.fmap_dims;
                beam.beamlet_size = args.beamlet_size;

                if (args.all_beamlets) {
                    beam.fluence_map = std::vector<float>(beam.fmap_size.x * beam.fmap_size.y, 1.0);
                } else {
                    beam.fluence_map = std::vector<float>(beam.fmap_size.x * beam.fmap_size.y);

                    // function implementation in "read_CT.cpp"
                    findFluenceProjection(
                            beam.fluence_map,
                            ptv_mask,
                            beam.isocenter,
                            beam.source,
                            beam.fmap_size,
                            beam.beamlet_size,
                            beam.azimuth,
                            beam.zenith,
                            beam.coll,
                            args.verbose
                            );

                    if (args.apertureready) {
                        fmap_post_apertureready(beam.fluence_map, beam.fmap_size);
                    }
                }

                // check beamlet counts for each beam
                int nbeamlets = 0;
                for (const& ff : beam.fluence_map) {
                    if (ff>0) {
                        nbeamlets++;
                        if (args.verbose<2) { break; }
                    }
                }
                if (args.verbose>=2) {
                    if (nbeamlets<=0) {
                        std::cout << set_color(COLOR::YELLOW);
                    }
                    std::cout << "  Beam..."<<std::setw(4)<<std::setfill('.')<<bb+1<<": #Beamlets: "<<nbeamlets << std::endl;
                    std::cout << set_color();
                } else if (nbeamlets<=0){
                    std::cout << set_color(COLOR::YELLOW) << "Beam #"<<bb+1<<" has 0 beamlets"<<set_color()<<std::endl;
                }
            }
        }
        logger.print_tail();
        std::cout << std::endl;

        if (args.verbose) {
            logger.print_head("FINAL BEAM SPECIFICATIONS");
            for (int b=0; b<beams.size(); b++) {
                std::cout << "  Beam..."<<std::setw(4)<<std::setfill('.')<<b+1<<": "<<beams[b] << std::endl;
            }
            logger.print_tail();
            std::cout << std::endl;
        }

        logger.print_head("LOADING MASK STRUCTURES");
        ROIMaskList roi_list {};

        // Fetch binary volumes for each ROI name
        for (const auto& roi_name : roi_names) {
            // Check for duplicates
            bool unique = true;
            for (const auto& r : roi_list.getROINames()) {
                if (r == roi_name) { unique = false; break; }
            }
            if (!unique) {
                std::cout << "Excluding redundant ROI specification: \""<<roi_name<< "\"" << std::endl<<std::endl;
                continue;
            }

            // validate name against rtstruct
            int roi_idx = getROIIndex(rtstruct, roi_name, true, args.verbose);
            StructureSet roi;
            if (!loadStructureSet(roi, rtstruct, roi_idx, args.verbose)) {
                if (!args.verbose) std::cout << "Failed to load ROI Data for: \""<< roi_name <<"\"" << std::endl;
                return 1;
            }

            // Construct BaseROIMask using rtstruct contour data from file
            ArrayProps roi_bbox = getROIExtents(roi, frameofref, args.verbose);
            Volume<uint8_t> roi_mask {};
            std::cout << "Creating ROI Mask" << std::endl;
            roi_mask = generateContourMask(roi, ctdata.get_frame(), density.get_frame());

            // crop mask to roi_bbox
            std::vector<uint8_t> cropped_mask(roi_bbox.nvoxels());
            for (uint ii=0; ii<roi_bbox.crop_size.x; ii++) {
                for (uint jj=0; jj<roi_bbox.crop_size.y; jj++) {
                    for (uint kk=0; kk<roi_bbox.crop_size.z; kk++) {
                        uint64_t full_key = ((kk+roi_bbox.crop_start.z)*roi_bbox.size.y + (jj+roi_bbox.crop_start.y))*roi_bbox.size.x + (ii+roi_bbox.crop_start.x); // iterator over roi_bbox volume
                        uint64_t crop_key = (kk*roi_bbox.crop_size.y + jj)*roi_bbox.crop_size.x + ii;
                        cropped_mask.at(crop_key) = roi_mask.at(full_key);
                    }
                }
            }

            // write binary volume to file to check
            if (debugwrite) {
                {
                    std::ostringstream outpath;
                    outpath << "mask_" << roi_name;
                    std::vector<float> temp(roi_mask.nvoxels());
                    for (uint i=0; i<roi_mask.nvoxels(); ++i) { temp[i] = roi_mask[i]; }
                    std::cout << set_color(COLOR::BLUE)<<"Writing mask to \""<<outpath.str()<<".raw"<<"\""<<set_color()<< std::endl;
                    write_debug_data<float>(temp.data(), roi_mask.size, outpath.str().c_str(), true);
                }
                {
                    std::ostringstream outpath;
                    outpath << "cropmask_" << roi_name;
                    std::vector<float> temp(cropped_mask.size());
                    for (uint i=0; i<cropped_mask.size(); ++i) { temp[i] = cropped_mask.at(i); }
                    std::cout << set_color(COLOR::BLUE)<<"Writing cropmask to \""<<outpath.str()<<".raw"<<"\""<<set_color()<< std::endl;
                    write_debug_data<float>(temp.data(), roi_bbox.crop_size, outpath.str().c_str(), true);
                }
            }

            roi_list.push_back(new DenseROIMask(roi_name, cropped_mask, roi_bbox));
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Discovered " << roi_list.size() << " ROIs"<< std::endl;
        uint idx = 0;
        for (const auto& v : roi_list.getROINames()) {
            ++idx;
            std::cout << "  " << idx << ": "<< v << std::endl;
        }
        std::ostringstream roi_list_outpath;
        roi_list_outpath << Paths::Instance()->temp_dir() << "/" << "roi_list.h5";
        if (args.verbose) { std::cout << "\nWriting ROI List to \""<<roi_list_outpath.str()<<"\"" << std::endl; }
        roi_list.writeToFile(roi_list_outpath.str());

        logger.print_tail();
        std::cout << std::endl;
    } else if (!args.densvol.empty()) {
        // unpack volume
        if (args.verbose) { printf("Reading Density from \"%s\"\n", args.densvol.c_str()); }
        FloatVolume vol;
        FloatVolume::readFromFile(vol, (const std::string)args.densvol);
        if (args.verbose) { std::cout << "Density:\n"<< vol << std::endl; }
        cudaCreateTexIso(density, vol, args.voxelsize, NULL, false);

        // use default bbox size
        calc_bbox_start = make_uint3(0,0,0);
        calc_bbox_size = density.size;
    } // else is handled in argparser

    std::cout << std::endl;
    logger.print_head("WRITING FLUENCE MAPS");
    uint bidx = 0;
    for (BEAM& beam : beams) {
        write_fluence_map( beam, bidx, beam.fmap_size, args.verbose>=2 );
        ++bidx;
    }
    logger.print_tail();
    std::cout << std::endl;

    if (write_omni_beam_list(beams, beams.size(), args.verbose) < 0) {
        printf("\n Failed to write omni beam list.\n");
        return -1;
    }

    if (args.timing) {
        printf("\nTime elapsed reading parameters & CT data: %4.3f msec\n", ((float)clock() - start_task)*1000 / CLOCKS_PER_SEC );
        start_task = clock();
    }

    //read spectrum data
    if ( (1 != read_spectrum_file(&mono_kernels,args.verbose)) ) {
        printf("Failed reading spectrum file!\n");
        exit(EXIT_FAILURE);
    }
    // monoenergetic kernel file
    if (args.verbose) {
        printf("Spectrum file is %s\n",mono_kernels.spectrum_file.c_str());
    }
        std::cout << poly_kernel.kernel_file<< std::endl;
    //read all monoenergetic kernels that were specified in spectrum file
    for (int i=0; i<mono_kernels.nkernels; i++) {
        if ( (1 != read_kernel(&(mono_kernels.kernel[i]))) ) {
            printf("Failed reading kernel!\n");
            exit(EXIT_FAILURE);
        }
    }
    //create polyenergetic kernel from mono kernels and fluence, mu data
    if ( (1 != make_poly_kernel(&mono_kernels,&poly_kernel,args.verbose,args.debug)) ) {
        printf("Failed making polyenergetic kernel!\n");
        exit(EXIT_FAILURE);
    }
    if (args.timing) {
        printf("\nTime elapsed during Read_Specfile, Read_Kernel, Make_Poly and Display_Info: %4.3f msec\n", ((float)clock() - start_task)*1000 / CLOCKS_PER_SEC );
        start_task = clock();
    }

    // write isotropic data volume out as a binary file
    write_debug_data<float>(density.data(), density.size, "density");
    // also write polyenergetic kernel and radial boundaries
    if ( (1 != make_cumulative_kernel(&poly_kernel, args.nphi, args.ntheta, args.verbose)) ) {
        printf("Failed creating cumulative poly-kernel\n");
        exit(EXIT_FAILURE);
    }

    // lastly write a header file with pertinent information on the data volume
    // such as size, resolution, nphi, ntheta, etc.
    {
        bool result = write_omni_header(
                density.size,
                density.voxsize,
                density.start,
                calc_bbox_start,
                calc_bbox_size,
                args.nphi, args.ntheta, N_KERNEL_RADII,
                beams.size(),
                args.beam_spec,
                args.ptv_name,
                args.rev_latspacing,
                args.rev_longspacing,
                args.kernel_extent,
                args.ss_factor,
                args.max_rev_size,
                args.penumbra
                );
        if (!result) {
            printf("Failed writing omni-header!\n");
            exit(EXIT_FAILURE);
        }

    }

    if (args.timing) {
        printf("\nTime elapsed during Data Write: %4.3f msec\n", ((float)clock() - start_task)*1000 / CLOCKS_PER_SEC );
        start_task = clock();
    }

    printf("\nTotal time elapsed: %4.3f msec\n\n", ((float)clock() - start_total) * 1000 / CLOCKS_PER_SEC);
    fcloseall();

    exit(EXIT_SUCCESS);
}
