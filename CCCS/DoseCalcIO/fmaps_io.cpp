#include "fmaps_io.h"

#include <iostream>
#include <list>
#include <helper_cuda.h>
#include <helper_math.h>
#include <memory>
#include <utility> // std::move

int _read_patient_metadata(H5::Group& h5group, HEADER_PATIENT& patient_header) {
    // read simple/scalar attributes
    {
        auto att = h5group.openAttribute("N_beams");
        att.read(H5::PredType::NATIVE_USHORT, &patient_header.N_beams);
    } {
        float temp[3];
        auto att = h5group.openAttribute("dicom_start_cm");
        att.read(H5::PredType::NATIVE_FLOAT, temp);
        ARR3VECT(patient_header.dicom_start_cm, temp);
    } {
        uint temp[3];
        auto att = h5group.openAttribute("full_dicom_size");
        att.read(H5::PredType::NATIVE_UINT, temp);
        ARR3VECT(patient_header.full_dicom_size, temp);
    } {
        float temp[3];
        auto att = h5group.openAttribute("voxel_size_cm");
        att.read(H5::PredType::NATIVE_FLOAT, temp);
        ARR3VECT(patient_header.voxel_size_cm, temp);
    } {
        uint temp[3];
        auto att = h5group.openAttribute("calc_bbox_start");
        att.read(H5::PredType::NATIVE_UINT, temp);
        ARR3VECT(patient_header.bbox_start, temp);
    } {
        uint temp[3];
        auto att = h5group.openAttribute("calc_bbox_size");
        att.read(H5::PredType::NATIVE_UINT, temp);
        ARR3VECT(patient_header.bbox_size, temp);
    } {
        auto att = h5group.openAttribute("sparsity_thresh");
        att.read(H5::PredType::NATIVE_FLOAT, &patient_header.sparsity_thresh);
    } {
        auto att = h5group.openAttribute("nphi");
        att.read(H5::PredType::NATIVE_USHORT, &patient_header.nphi);
    } {
        auto att = h5group.openAttribute("ntheta");
        att.read(H5::PredType::NATIVE_USHORT, &patient_header.ntheta);
    } {
        auto att = h5group.openAttribute("nradii");
        att.read(H5::PredType::NATIVE_USHORT, &patient_header.nradii);
    } {
        auto att = h5group.openAttribute("convlat_cm");
        att.read(H5::PredType::NATIVE_FLOAT, &patient_header.rev_latspacing_cm);
    } {
        auto att = h5group.openAttribute("convstep_cm");
        att.read(H5::PredType::NATIVE_FLOAT, &patient_header.rev_longspacing_cm);
    } {
        auto att = h5group.openAttribute("beam_spectrum");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        patient_header.beam_spectrum = buf;
    } {
        auto att = h5group.openAttribute("target_structure");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        patient_header.target_structure = buf;
    }
    // only for reduced case
    if (h5group.attrExists("roi_order")) {
        auto att = h5group.openAttribute("roi_order");
        hsize_t N = att.getSpace().getSimpleExtentNpoints();
        std::unique_ptr<char*[]> temp(new char*[N]);
        att.read(att.getDataType(), (void*)temp.get());
        patient_header.roi_order = std::vector<std::string>(&temp[0], &temp[N]);
    }
    if (h5group.attrExists("row_block_capacities")) {
        auto att = h5group.openAttribute("row_block_capacities");
        hsize_t N = att.getSpace().getSimpleExtentNpoints();
        std::unique_ptr<uint64_t[]> temp(new uint64_t[N]);
        att.read(H5::PredType::NATIVE_UINT64, temp.get());
        patient_header.row_block_capacities = std::vector<uint64_t>(&temp[0], &temp[N]);
    }

    return true;
}
int _read_beams(H5::Group& h5group, std::vector<HEADER_BEAM>& beam_headers) {
    // read beams, their metadata, and fluence maps from h5 file
    struct IndexedBeamHeader {
        IndexedBeamHeader(uint16_t idx, HEADER_BEAM&& beam_header) : idx{idx} { this->beam_header = std::move(beam_header); }
        uint16_t idx;
        HEADER_BEAM beam_header;
    };
    struct OpData { OpData(H5::Group& g) : h5group{g} {}; std::list<IndexedBeamHeader> groups={}; H5::Group& h5group; };

    OpData opdata{h5group};
    int iter_idx = 0; // iter_count is returned here
    h5group.iterateElems(".", &iter_idx,
        [](hid_t loc_id, const char* name, void* opdata) -> herr_t {
            // iterator body
            OpData* data = static_cast<OpData*>(opdata);
            H5::Group beam_group = data->h5group.openGroup(name);
            HEADER_BEAM this_beam_header {};
            // read beam compound datatype
            BEAM::_readFromHDF5(this_beam_header.beam_specs, beam_group);

            // read beam_header attributes
            { // get uid
                auto att = beam_group.openAttribute("beam_uid");
                att.read(H5::PredType::NATIVE_UINT, &this_beam_header.beam_uid);
                this_beam_header.beam_specs.uid = this_beam_header.beam_uid;
            } {
                auto att = beam_group.openAttribute("N_beamlets");
                att.read(H5::PredType::NATIVE_USHORT, &this_beam_header.N_beamlets);
            }
            data->groups.push_back( IndexedBeamHeader(this_beam_header.beam_uid, std::move(this_beam_header)) );

            return 0;
        }, (void*)&opdata);

    // sort headers by uid
    opdata.groups.sort( [](IndexedBeamHeader& a, IndexedBeamHeader& b) -> bool {
            // true if a belongs before b
            return (a.idx <= b.idx);
            } );

    // move headers to output vector
    for (auto&& indexed_beam_header : opdata.groups) {
        beam_headers.push_back( indexed_beam_header.beam_header ); // move, rvalue-ref
    }

    return true;
}

int read_fmaps(const std::string& filename, HEADER_PATIENT* patient_header, std::vector<HEADER_BEAM>& beam_headers ) {
    // read output of fluence map optimizer and return structs containing calculation metadata, beam specs,
    // and fluence intensities
    H5::H5File file(filename, H5F_ACC_RDONLY);

    // read patient metadata
    if (patient_header != nullptr) {
        H5::Group ptmeta_grp = file.openGroup("calc_specs");
        if (!_read_patient_metadata(ptmeta_grp, *patient_header)) { return false; }
    }

    // read beams
    H5::Group beams_grp = file.openGroup("beams");
    // test for new format
    try {
        beams_grp = beams_grp.openGroup("metadata");
    } catch (H5::GroupIException) { }
    if (!_read_beams(beams_grp, beam_headers)) { return false; }

    return true;
}
int read_fmaps(const std::string& filename, HEADER_PATIENT* patient_header, std::vector<BEAM>& beams) {
    std::vector<HEADER_BEAM> beam_headers;
    if (!read_fmaps(filename, patient_header, beam_headers)) { return false; }
    // extract beams vector
    for (auto&& beam_header : beam_headers) {
        beams.push_back( beam_header.beam_specs );
    }
    return true;
}

int _write_file_version(H5::Group& h5group, uint ftmagic, uint ftversionmajor, uint ftversionminor) {
    auto fgroup = h5group.createGroup("filetype");
    H5::DataSpace scalarspace;
    {
        auto att = fgroup.createAttribute("ftmagic", H5::PredType::STD_U8LE, scalarspace);
        att.write(H5::PredType::NATIVE_UINT, &ftmagic);
    }
    {
        auto att = fgroup.createAttribute("ftversionmajor", H5::PredType::STD_U8LE, scalarspace);
        att.write(H5::PredType::NATIVE_UINT, &ftversionmajor);
    }
    {
        auto att = fgroup.createAttribute("ftversionminor", H5::PredType::STD_U8LE, scalarspace);
        att.write(H5::PredType::NATIVE_UINT, &ftversionminor);
    }
    return 1;
}
