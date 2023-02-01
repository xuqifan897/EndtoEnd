#include "sparsify_worker.h"

#include <iostream>
#include <sstream>

SparsifyWorkerBase::~SparsifyWorkerBase() { s_workerid--; }
void SparsifyWorkerBase::join() {
    // defined to ensure deactivation of wqueue to kill wqueue after emptied before join takes place
    m_queue.deactivate();
    Thread::join();
}



void BeamletSparsifyWorker::run() {
    SparsifyData data;
    while(true) {
        if (m_queue.try_pop(data, workerid())) { // received valid data
            WLOG(m_workerid, "Preparing to write to file");
            WLOG(m_workerid, "recieved data: " << &data);
            write_sparse_beamlet_to_file(data.filename, data.dense_array, data.props, sparse_to_hdf5, (void*)&data.beamlet_header);
            WLOG(m_workerid, "done writing to file");
            if (data.dense_array != nullptr) {
                delete [] data.dense_array;
                data.dense_array = nullptr;
            }
        } else { WLOG(m_workerid, "Killing worker"); return; /* recieved kill signal */ }
    }
}


FullSparsifyWorker::FullSparsifyWorker(wqueue<SparsifyData>& queue, const std::string h5filename, HEADER_PATIENT* patient_header, const float sparsity_thresh) : SparsifyWorkerBase(queue, sparsity_thresh), m_h5filename{h5filename}, patient_header{patient_header} {
    // open central hdf5 file to use for life of the worker
    m_h5file = H5::H5File(m_h5filename, H5F_ACC_TRUNC);

    // Create secondary wqueue just for writing to file,
}
int FullSparsifyWorker::write_sparse_full(H5::H5File& h5file, SparsifyData& sdata) {
    // Write all beamlets added to wqueue to a single hdf5 file organized in hierarchy of beam_uid -> beamlet_uid -> sparse coeff data
    SparseData data;
    _sparsify_and_reduce(data, sdata.dense_array, sdata.props, sdata.roi_list, m_thresh);

    // store discovered data in header
    HEADER_BEAMLET& beamlet_header = sdata.beamlet_header;
    beamlet_header.N_coeffs = data.size();
    beamlet_header.row_block_sizes = data.row_block_sizes;

    // open beam group
    ushort beam_uid = sdata.beam_header.beam_uid;
    std::ostringstream groupname;
    groupname << "/beam_" << beam_uid;

    H5::Group beamgroup;
    if (m_beams_seen.find(beam_uid) == m_beams_seen.end()) {
        m_beams_seen.emplace(beam_uid);
        WLOG(m_workerid, "First time seeing beam: " << beam_uid);
        // This is first time beam has been seen
        WLOG(m_workerid, "creating group: " << "\"" << groupname.str() << "\"");
        beamgroup = h5file.createGroup(groupname.str());
        WLOG(m_workerid, "done creating group");
        // also write beam metadata
        WLOG(m_workerid, "writing beam header");
        _write_beam_metadata(beamgroup, sdata.beam_header);
        WLOG(m_workerid, "done writing beam header");
    } else {
        WLOG(m_workerid, "already seen beam: " << beam_uid);
        WLOG(m_workerid, "opening group: " << "\"" << groupname.str() << "\"");
        beamgroup = h5file.openGroup(groupname.str());
        WLOG(m_workerid, "done opening group");
    }

    // create beamlet group
    groupname << "/" << "beamlet_" << sdata.beamlet_header.beamlet_uid;
    WLOG(m_workerid, "creating group: " << "\"" << groupname.str() << "\"");
    auto beamletgroup = h5file.createGroup(groupname.str());
    WLOG(m_workerid, "done creating group");
    WLOG(m_workerid, "writing beamlet to file");
    _write_beamlet_to_hdf5(beamletgroup, data, sdata.beamlet_header, get_compress_lvl(), get_chunksize());
    WLOG(m_workerid, "done writing beamlet to file");

    return 1;
}
void FullSparsifyWorker::run() {
    // write patient metadata to root group once
    H5::Group rootgroup = m_h5file.openGroup("/");
    _write_patient_metadata(rootgroup, *patient_header);

    while(true) {
        SparsifyData data {};
        if (m_queue.try_pop(data, workerid())) { // received valid data
            WLOG(m_workerid, "Preparing to write to file");
            WLOG(m_workerid, "recieved data: " << &data);
            write_sparse_full(m_h5file, data);
            WLOG(m_workerid, "done writing to file");
            if (data.dense_array != nullptr) {
                WLOG(m_workerid, "Deleting " << data.dense_array);
                delete [] data.dense_array;
                data.dense_array = nullptr;
            } else {
                WLOG(m_workerid, "Not deleting " << data.dense_array);
            }
        } else { WLOG(m_workerid, "Killing worker"); return; /* recieved kill signal */ }
    }
}

// zero-init static variable
unsigned int SparsifyWorkerBase::s_workerid = 0;
