#include "sparse_io.h"

#include <cstdio>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <map>
#include <algorithm>
#include <exception>
#include <cassert>
#include <memory>

#include "version.h"
#include "dosecalc_defs.h"
#include "./io_helpers.h"
#include "./macros.h"
#include "Utilities/timing.h"

using namespace dcio;

void _sparsify(SparseData& sparsedata, const float* array, ArrayProps& props, float thresh) {
    if (props.crop_size.x <= 0 || props.crop_size.y <= 0 || props.crop_size.z <= 0) {
        props.crop_size = props.size;
    }

    // density statistics
    uint64_t nnonzero = 0;
    uint64_t ndropped = 0;
    uint64_t ntotal = props.nvoxels();

    float max_val = *std::max_element(array, array+props.nvoxels());
    // array looping - sparsification
    for (uint i=0; i<props.crop_size.x; i++) {
        for (uint j=0; j<props.crop_size.y; j++) {
            for (uint k=0; k<props.crop_size.z; k++) {
                uint64_t idx = (k*props.crop_size.y + j)*props.crop_size.x + i;

                // eval non-zero
                float val = array[idx];
                // std::cout << "idx: " << idx << " - val: " << std::setprecision(3) << val << std::endl;
                if (val <= thresh*max_val) {
                    if (val>0) { ndropped++; }
                    continue;
                }

                // store kvp to vector (implicitly sorted by looping order)
                uint64_t true_idx = ((k+props.crop_start.z)*props.size.y + (j+props.crop_start.y))*props.size.x + (i+props.crop_start.x);
                // vect.push_back( KeyValPair{true_idx, val} );
                sparsedata.kvp.keys.push_back(true_idx);
                sparsedata.kvp.vals.push_back(val);
                nnonzero++;
            }
        }
    }

    sparsedata.perc_nonzero = 100.0*nnonzero/ntotal;
    sparsedata.perc_dropped = 100.0*ndropped/ntotal;
}

//TODO: pre-cache mask indices rather than recompute for each beamlet
void _reduce_by_roi(SparseData& reduceddata, SparseData& sparsedata, ROIMaskList& roi_list, ArrayProps& props) {
     /* Accepts sparsedata as a sparse column of keys (over full dicom volume) and associated dose contribs.
      * and produces reduceddata which is sparsedata with rows removed where the corresponding voxel is not
      * present in any of the input ROIs and thus has no effect on the optimization result
      *
      * reduceddata is further ordered such that PTV is the first row-block, followed by each OAR in arbitrary
      * but defined order. The key for each block becomes the index into the ROI not including any voxels that
      * are masked out of the ROI. the first possible key value for each block depends on the total number of
      * potential rows in the previous blocks (sum of total unmasked voxels in all previous ROIs)
      *
      * outputs:
      *     - reduceddata [KeyValPairs]: key is converted to encode which ROI and index of voxel within that
      *         ROI to which the accompanying dose-value corresponds.
      *     - ordered list of ROIs and accompanying Nvoxel counts. (sum of Nvoxel counts should be total #rows
      *         in dense representation of sparse column vector).
      * */

    // create sparsevect for efficient random access based on key (time: 3ms, mem: 64-500MB)
    // AutoTimer timer_task;
    std::vector<float> sparsevect(props.size.x*props.size.y*props.size.z, 0.f);
    for (uint64_t i=0; i<sparsedata.kvp.keys.size(); ++i) {
        uint64_t key = sparsedata.kvp.keys.at(i);
        sparsevect.at(key) = sparsedata.kvp.vals.at(i);
    }
    // timer_task.stop_print_time_elapsed("sparsevect creation");

    // Guess at appropriate size for reduceddata keys/vals (time: 0ms)
    uint64_t MAGIC = static_cast<uint64_t>(
            sparsedata.size() * 0.4 +                            // PTV nonzero Density
            (roi_list._coll.size()-1)*sparsedata.size() * 0.15   // OAR nonzero Density
            );
    reduceddata.kvp.keys.reserve( MAGIC );
    reduceddata.kvp.vals.reserve( MAGIC );
    // printf("roi_list.size = %zu\n", roi_list._coll.size());
    // printf("sparsedata.size = %zu\n", sparsedata.size());
    // printf("MAGIC = %lu\n", MAGIC);
    // printf("capacity: %lu\n", reduceddata.keys.capacity());

    // density statistics
    uint64_t nnonzero = 0;
    uint64_t ntotal = props.nvoxels();

    // Loop over ROI structures within calc_bbox
    uint64_t mark = 0;
    for (auto &roi : roi_list._coll) {
        // ensure mask dims match Dicom dims (so calc_bbox spec can be applied to both)
        if (!(props.size.x == roi->props.size.x && props.size.y == roi->props.size.y && props.size.z == roi->props.size.z)) {
            throw std::runtime_error("props.size != roi.props.size");
        }
        // printf("roi: %s, nones: %d\n",roi->name.c_str(),roi->nones());

        // iterate over array and keep elements where mask == 1
        uint64_t Nkept = 0;
        uint64_t elem_key = 0; // iterator only over voxels in ROI
        // Implement bbox over just the ROI volume instead of full dose calc volume
        // This ensures that every voxel in ROI is touched and that iterator has fewer points to test
        for (uint ii=0; ii<roi->props.crop_size.x; ii++) {
            for (uint jj=0; jj<roi->props.crop_size.y; jj++) {
                for (uint kk=0; kk<roi->props.crop_size.z; kk++) {
                    uint64_t crop_key = (kk*roi->props.crop_size.y + jj)*roi->props.crop_size.x + ii;
                    uint64_t orig_key = ((kk+roi->props.crop_start.z)*roi->props.size.y + (jj+roi->props.crop_start.y))*roi->props.size.x + (ii+roi->props.crop_start.x); // iterator over roi_bbox volume

                    // convert to reduced matrix key
                    // printf("crop_key: %lu, bool: %u\n", crop_key, (*roi)[crop_key]);
                    if ((*roi)[crop_key]) {
                        float& dose_val = sparsevect.at(orig_key);
                        if (dose_val) { // || sparsemap.count(orig_key)) {
                            uint64_t pos = mark + elem_key;
                            // if (pos>=max_key) {
                            //     // std::cout << set_color(COLOR::RED) << "ASSERT: max_key:"<<max_key<<" roi: "<<roi->name<< " ii:"<<ii<<" jj:"<<jj<<" kk:"<<kk<<" crop_key:"<<crop_key<<" orig_key:"<<orig_key<<" crop_size:("<<props.crop_size.x<<","<<props.crop_size.y<<","<<props.crop_size.z<<") crop_start:("<<props.crop_start.x<<","<<props.crop_start.y<<","<<props.crop_start.z<<")"<<" size:("<<props.size.x<<","<<props.size.y<<","<<props.size.z<<" elem_key:"<<elem_key<<" mark:"<<mark<<" pos:"<<pos<<" dose_val:"<<dose_val<<" nkept:"<<Nkept<<set_color()<<std::endl;
                            //     throw std::runtime_error("stored sparse key was greater than expected");
                            // }
                            // printf("%lu, %0.2e\n", pos, sparsevect.at(orig_key));
                            reduceddata.kvp.keys.push_back(pos);
                            reduceddata.kvp.vals.push_back(dose_val);
                            ++Nkept;
                        }
                        ++elem_key;
                    }
                }
            }
        }
        // printf("#elements in mask(%d) || #elements in reduced vector for ROI=\"%s\" (%d)\n", roi->nones(), roi->name.c_str(), Nkept);
            if (elem_key != roi->nones()) {
                throw std::runtime_error("elem_key != roi->nones()");
            }
        assert(elem_key == roi->nones() /* "# iterated roi elements must equal number of masked ROI voxels" */);
        mark += elem_key;
        reduceddata.row_block_sizes.push_back(Nkept);
        nnonzero += Nkept;
    }
    double perc_nonzero = 100.0*nnonzero/ntotal;
    reduceddata.perc_reducedrop = sparsedata.perc_nonzero-perc_nonzero;
    reduceddata.perc_nonzero = perc_nonzero;
    reduceddata.perc_dropped = sparsedata.perc_dropped;
    // printf("total elements in reduced: %lu\n", reduceddata.keys.size());
}
void _sparsify_and_reduce(SparseData& data, const float* array, ArrayProps& props, ROIMaskList* roi_list, float thresh) {
    // ALWAYS CALL THIS for sparsify task. If roi_list is non-empty, reduction and reordering of rows will
    // also occur
    if (roi_list != nullptr && roi_list->size()) {
        // store to tempvar for use in reduce
        SparseData sparsedata;
        // AutoTimer timer_task;
        _sparsify(sparsedata, array, props, thresh);
        // timer_task.restart_print_time_elapsed("sparsify");
        _reduce_by_roi(data, sparsedata, *((ROIMaskList*)roi_list), props);
        // timer_task.stop_print_time_elapsed("reduce");
    } else {
        // directly modify outvar
        // AutoTimer timer_task;
        _sparsify(data, array, props, thresh);
        // timer_task.stop_print_time_elapsed("sparsify");
    }

}
int _write_beamlet_to_hdf5(H5::Group& h5group, const SparseData& data, HEADER_BEAMLET& beamlet_header, int compress_lvl, uint64_t chunksize) {
    hsize_t dims[] = { beamlet_header.N_coeffs };
    H5::DataSpace simplespace(1, dims);
    H5::DataSpace scalarspace;

    H5::DSetCreatPropList plist = static_cast<H5::DSetCreatPropList>(H5::DSetCreatPropList::DEFAULT);
    if (compress_lvl >= 0) {
        plist = H5::DSetCreatPropList();
        hsize_t chunk_dims[] = { std::min(beamlet_header.N_coeffs, chunksize) };
        plist.setChunk(1, chunk_dims);
        plist.setDeflate(std::min(9, compress_lvl));
    }

    // Create dataset - keys
    {
        std::ostringstream dname;
        dname << "lindex";
        auto dset_keys = h5group.createDataSet(dname.str().c_str(), H5::PredType::STD_U64LE, simplespace, plist);
        dset_keys.write((void*)&data.kvp.keys.front(), H5::PredType::NATIVE_UINT64);
    }

    // Create dataset - values
    {
        std::ostringstream dname;
        dname << "coeffs";
        auto dset_vals = h5group.createDataSet(dname.str().c_str(), H5::PredType::IEEE_F32LE, simplespace, plist);
        dset_vals.write((void*)&data.kvp.vals.front(), H5::PredType::NATIVE_FLOAT);
    }

    return 1;
}
int _write_beamlet_metadata(H5::Group& h5group, const HEADER_BEAMLET& beamlet_header, const SparseData& data) {
    // write attributes
    H5::DataSpace scalarspace;
    {
        auto att = h5group.createAttribute("beamlet_uid", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &beamlet_header.beamlet_uid);
    } {
        auto att = h5group.createAttribute("N_coeffs", H5::PredType::STD_U64LE, scalarspace);
        att.write(H5::PredType::NATIVE_UINT64, &beamlet_header.N_coeffs);
    } {
        auto att = h5group.createAttribute("perc_nonzero", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &data.perc_nonzero);
    } {
        auto att = h5group.createAttribute("perc_dropped", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &data.perc_dropped);
    }
    // only for reduced case
    if (beamlet_header.reduced()) {
        // write attribute - ord. list of # nonzeros per row-block
        hsize_t dims[] = {beamlet_header.row_block_sizes.size()};
        auto att = h5group.createAttribute("row_block_sizes", H5::PredType::STD_U64LE, H5::DataSpace(1, dims));
        std::vector<uint64_t> temp(beamlet_header.row_block_sizes.size());
        uint ii = 0;
        for (const auto& xx : beamlet_header.row_block_sizes) {
            temp[ii] = xx;
            ++ii;
        }
        att.write(H5::PredType::NATIVE_UINT64, temp.data());
        {
            auto att = h5group.createAttribute("perc_reducedrop", H5::PredType::IEEE_F32LE, scalarspace);
            att.write(H5::PredType::NATIVE_FLOAT, &data.perc_reducedrop);
        }
    }

    return 1;
}
int _write_beam_metadata(H5::Group& h5group, const HEADER_BEAM& beam_header) {
    // Write beam header to hdf5 group as attributes
    H5::DataSpace scalarspace;

    // write simple/scalar attributes
    {
        auto att = h5group.createAttribute("beam_uid", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &beam_header.beam_uid);
    } {
        auto att = h5group.createAttribute("N_beamlets", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &beam_header.N_beamlets);
    }

    // write beam data as compound datatype
    beam_header.beam_specs._writeToHDF5(h5group);

    return 1;
}
int _write_patient_metadata(H5::Group& h5group, const HEADER_PATIENT& patient_header) {
    h5group = h5group.createGroup("calc_specs");

    H5::DataSpace scalarspace;

    // tuple dataspaces
    hsize_t tuple3_dims[] = { 3 };

    // write simple/scalar attributes
    {
        std::string vers = VERSION_STRING;
        H5::StrType str_t{H5::PredType::C_S1, vers.length()+1};
        auto att = h5group.createAttribute("dosecalc_version", str_t, scalarspace);
        att.write(str_t, vers);
    } {
        auto att = h5group.createAttribute("N_beams", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &patient_header.N_beams);
    } {
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_FLOAT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::IEEE_F32LE, 1, tuple3_dims);
        float temp[3];
        VECT3ARR(temp, patient_header.dicom_start_cm);
        auto att = h5group.createAttribute("dicom_start_cm", tuple3_t, scalarspace);
        att.write(tuple3_native_t, temp);
    } {
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_UINT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::STD_U16LE, 1, tuple3_dims);
        uint temp[3];
        VECT3ARR(temp, patient_header.full_dicom_size);
        auto att = h5group.createAttribute("full_dicom_size", tuple3_t, scalarspace);
        att.write(tuple3_native_t, temp);
    } {
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_FLOAT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::IEEE_F32LE, 1, tuple3_dims);
        float temp[3];
        VECT3ARR(temp, patient_header.voxel_size_cm);
        auto att = h5group.createAttribute("voxel_size_cm", tuple3_t, scalarspace);
        att.write(tuple3_native_t, temp);
    } {
        auto att = h5group.createAttribute("convlat_cm", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &patient_header.rev_latspacing_cm);
    } {
        auto att = h5group.createAttribute("convstep_cm", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &patient_header.rev_longspacing_cm);
    } {
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_UINT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::STD_U16LE, 1, tuple3_dims);
        uint temp[3];
        VECT3ARR(temp, patient_header.bbox_start);
        auto att = h5group.createAttribute("calc_bbox_start", tuple3_t, scalarspace);
        att.write(tuple3_native_t, temp);
    } {
        H5::ArrayType tuple3_native_t(H5::PredType::NATIVE_UINT, 1, tuple3_dims);
        H5::ArrayType tuple3_t(H5::PredType::STD_U16LE, 1, tuple3_dims);
        uint temp[3];
        VECT3ARR(temp, patient_header.bbox_size);
        auto att = h5group.createAttribute("calc_bbox_size", tuple3_t, scalarspace);
        att.write(tuple3_native_t, temp);
    } {
        auto att = h5group.createAttribute("penumbra_cm", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &patient_header.penumbra_cm);
    } {
        auto att = h5group.createAttribute("sparsity_thresh", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &patient_header.sparsity_thresh);
    } {
        auto att = h5group.createAttribute("nphi", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &patient_header.nphi);
    } {
        auto att = h5group.createAttribute("ntheta", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &patient_header.ntheta);
    } {
        auto att = h5group.createAttribute("nradii", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_USHORT, &patient_header.nradii);
    } {
        auto att = h5group.createAttribute("kernel_extent_cm", H5::PredType::IEEE_F32LE, scalarspace);
        att.write(H5::PredType::NATIVE_FLOAT, &patient_header.kernel_extent_cm);
    } {
        H5::StrType str_t{H5::PredType::C_S1, patient_header.beam_spectrum.length()+1};
        auto att = h5group.createAttribute("beam_spectrum", str_t, scalarspace);
        att.write(str_t, patient_header.beam_spectrum);
    } {
        H5::StrType str_t{H5::PredType::C_S1, patient_header.target_structure.length()+1};
        auto att = h5group.createAttribute("target_structure", str_t, scalarspace);
        att.write(str_t, patient_header.target_structure);
    }
    // only for reduced case
    if (patient_header.reduced()) {
        // write attribute - ord. list of capacities per block
        {
            hsize_t dims[] = {patient_header.row_block_capacities.size()};
            auto att = h5group.createAttribute("row_block_capacities", H5::PredType::STD_U64LE, H5::DataSpace(1, dims));
            att.write(H5::PredType::NATIVE_UINT64, patient_header.row_block_capacities.data());
        }
        // write attribute - ord. list of roi name strings
        {
            hsize_t dims[] = {patient_header.roi_order.size()};
            H5::StrType str_t(H5::PredType::C_S1, H5T_VARIABLE);
            auto att = h5group.createAttribute("roi_order", str_t, H5::DataSpace(1, dims));
            std::vector<const char*> temp;
            for (auto& xx : patient_header.roi_order) {
                temp.push_back( xx.c_str() );
            }
            att.write(str_t, temp.data());
        }
    }
    return 1;
}

int sparse_to_text(const std::string& filename, const SparseData& data, const void* _) {
    // defines method to store all data from sparse vector to text file

    // add extension
    std::ostringstream fullfilename;
    fullfilename << filename << ".txt";

    std::ofstream outfile(fullfilename.str());
    outfile.setf(std::ios::scientific, std::ios::floatfield);
    int prec = 6;
    for (uint iii=0; iii<data.kvp.keys.size(); iii++) {
        outfile << std::setprecision(0) << data.kvp.keys[iii] << " " << std::setprecision(prec) << data.kvp.vals[iii] << std::endl;
    }
    outfile.close();

    return 1;
}
int sparse_to_binary(const std::string& filename, const SparseData& data, const void* _) {
    // defines method to store all data from sparse vector to binary file
    // Data stored as raw bytes with no metadata or header - alternating between:
    // -- 4 byte (uint):  linearized index into full dicom array
    // -- 4 byte (float): dose coeff. value; edge between beamlet and dicom voxel


    // add extension
    std::ostringstream fullfilename;
    fullfilename << filename << ".bin";

    FILE* outfile;
    outfile = fopen(fullfilename.str().c_str(), "wb");
    fwrite((char*)&data.kvp.keys[0], sizeof(uint), data.kvp.keys.size(), outfile);
    fwrite((char*)&data.kvp.vals[0], sizeof(float), data.kvp.vals.size(), outfile);
    fclose(outfile);

    return 1;
}
int sparse_to_hdf5(const std::string& filename, const SparseData& data, const void* header) {
    // add extension
    std::ostringstream fullfilename;
    fullfilename << filename << ".h5";

    // Create file
    H5::H5File file(fullfilename.str(), H5F_ACC_TRUNC);
    H5::Group rootgroup = file.openGroup("/");

    //write
    _write_beamlet_to_hdf5(rootgroup, data, *((HEADER_BEAMLET*)header));

    return 1;
}


int write_sparse_beamlet_to_file(const std::string& filename, const float* array, ArrayProps& props, int(*tofile)(const std::string&, const SparseData&, const void* header), void* header, ROIMaskList* roi_list, const float thresh) {
    // loop through dense array of size: arr_size, encoding non-zero vals in Dict of Keys style sparse matrix
    SparseData data;
    _sparsify_and_reduce(data, array, props, roi_list, thresh);

    // store discovered data in header
    if (header!=nullptr) {
        HEADER_BEAMLET* beamlet_header = (HEADER_BEAMLET*)header;
        beamlet_header->N_coeffs = data.size();
        beamlet_header->row_block_sizes = data.row_block_sizes;
    }

    // write to file
    if (!tofile(filename, data, header)) {
        std::cout << "Failed storing sparse data to file: \""<< filename <<"\"" << std::endl;
        return 0;
    }

    return 1; //success
}
int write_beam_metadata(const std::string& filename, const HEADER_BEAM& beam_header) {
    // Write beam header data to file in sparse data directory
    // add extension
    std::ostringstream fullfilename;
    fullfilename << filename << ".h5";

    H5::H5File file(fullfilename.str(), H5F_ACC_TRUNC);
    H5::Group rootgroup = file.openGroup("/");
    _write_beam_metadata(rootgroup, beam_header);

    return 1;
}
int write_patient_metadata(const std::string& filename, const HEADER_PATIENT& patient_header) {
    // Write patient header data to file in sparse data directory

    // add extension
    std::ostringstream fullfilename;
    fullfilename << filename << ".h5";

    H5::H5File file(fullfilename.str(), H5F_ACC_TRUNC);
    H5::Group rootgroup = file.openGroup("/");
    _write_patient_metadata(rootgroup, patient_header);

    return 1;
}


int test_read_binary(const std::string& filename, KeyValPairs& in_sparse) {
    uint64_t fsize = filesize(filename);
    // std::cout << "Filesize: " << fsize << " bytes" << std::endl;

    std::ifstream infile(filename, std::ios::binary);
    char* fBuf = new char[fsize];
    infile.read((char*)fBuf, fsize);

    uint64_t N = fsize / (sizeof(uint64_t)+sizeof(float));

    // std::cout << "#elements: " << N << " bytes" << std::endl;

    in_sparse.keys.clear(); in_sparse.vals.clear();
    for (uint64_t ii=0; ii<N; ii++) {
        in_sparse.keys.push_back( *((uint64_t*)(fBuf + ii*sizeof(uint64_t))) );
    }
    uint64_t offset = N*sizeof(uint64_t);
    for (uint64_t ii=0; ii<N; ii++) {
        in_sparse.vals.push_back( *((float*)(fBuf + offset + ii*sizeof(float))) );
    }
    // std::cout << "Reading in sparse vector of size: " << in_sparse.size() << std::endl;

    infile.close();

    return 1;

}
