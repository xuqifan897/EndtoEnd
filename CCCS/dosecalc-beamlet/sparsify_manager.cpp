#include "sparsify_manager.h"

#include <exception>
#include <utility>
#include "DoseCalcIO/fmaps_io.h"
#include "DoseCalcIO/io_helpers.h"
#include "Utilities/timing.h"
#include "./debug.h"
#include "hdf5_hl.h"

static uint FTMAGIC = 0x2B;
static uint FTVERSIONMAJOR = 1;
static uint FTVERSIONMINOR = 2;

void SRWorker::run() {
    while(true) {
        SRData data {}; // data must go out of scope each loop for dense_array to be destroyed
        if (parent.m_srqueue.try_pop(data, m_workerid)) { // recvd valid data via move assignment
            WLOG(m_workerid, "Preparing to sparsify+reduce");
            WLOG(m_workerid, "recieved data: " << &data);
            SparseData sparsedata {};
            if (data.dense_array) {
                _sparsify_and_reduce(sparsedata, data.dense_array.get(), data.props, data.roi_list, parent.m_thresh);
                WLOG(m_workerid, "done sparsify+reduce");
            }

            // package and send to w2fqueue via move
            WLOG(m_workerid, "Preparing to push sparsedata");
            BaseW2FWorker::W2FData w2fdata{};
            w2fdata.beam_header = std::move(data.beam_header);
            if (sparsedata.size()) {
                w2fdata.sparsedata = std::move(sparsedata);
                w2fdata.beamlet_header = std::move(data.beamlet_header);
            }
            parent.m_w2fqueue.push(std::move(w2fdata));
            WLOG(m_workerid, "Done pushing sparsedata");

            // return dense_array pointer to the stack
            WLOG(m_workerid, "returning dense_array block to stack");
            if (data.dense_array) {
                memset(data.dense_array.get(), 0, data.props.nvoxels());
                parent.push_memblock(data.dense_array);
            }
        } else { WLOG(m_workerid, "Killing worker"); return; /* recvd kill-signal from wqueue */ }
    }
}
void BaseW2FWorker::run() {
    while (true) {
        BaseW2FWorker::W2FData data {}; // data must go out of scope each loop for sparsedata to be destroyed
        if (parent.m_w2fqueue.try_pop(data)) { // recvd valid data via move assign
            // store discovered meta in header
            data.beamlet_header.N_coeffs = data.sparsedata.size();
            data.beamlet_header.row_block_sizes = data.sparsedata.row_block_sizes;

            // write_to_file(data);
            write_beam_metadata(data);
            if (data.beamlet_header.N_coeffs>0) {
                write_beamlet(data);
            }
        } else { WLOG(m_workerid, "Killing worker"); return; /* recvd kill-signal from wqueue */ }
    }
};
void BaseW2FWorker::write_beam_metadata(BaseW2FWorker::W2FData& data) {
    unsigned short& beam_uid = data.beam_header.beam_uid;
    std::unique_lock<std::mutex> lock(s_central_mutex);
    if (s_beams_seen.count(beam_uid) == 0) {
        s_beams_seen.emplace(data.beam_header.beam_uid);
        WLOG(m_workerid, "First time seeing beam: " << beam_uid);

        // open beam group
        std::ostringstream datagroupname;
        datagroupname << SparsifyManager::beams_data_group() + "/beam_" << std::setfill('0') << std::setw(5) << beam_uid;
        parent.m_h5file.createGroup(datagroupname.str());
        WLOG(m_workerid, "creating group: " << "\"" << datagroupname.str() << "\"");

        std::ostringstream metagroupname;
        metagroupname << SparsifyManager::beams_meta_group() + "/beam_" << std::setfill('0') << std::setw(5) << beam_uid;
        WLOG(m_workerid, "creating group: " << "\"" << metagroupname.str() << "\"");
        H5::Group beamgroup = parent.m_h5file.createGroup(metagroupname.str());
        WLOG(m_workerid, "done creating group");

        // also write beam metadata
        WLOG(m_workerid, "writing beam header");
        _write_beam_metadata(beamgroup, data.beam_header);
        WLOG(m_workerid, "done writing beam header");

        // create per-beam output sidecar directory
        if (parent.m_write_strategy == WRITE_STRATEGY::PER_BEAMLET) {
            std::ostringstream beam_sidecar_path;
            beam_sidecar_path << parent.m_sidecar_path << "beam_" << std::setfill('0') << std::setw(5) << data.beam_header.beam_uid << std::setfill('\0') << "/";
            WLOG(m_workerid, "Creating directory: " << beam_sidecar_path.str());
            // if (dcio::dir_exists(beam_sidecar_path.str())) {
            //     dcio::remove_directory(beam_sidecar_path.str(), false, true);
            // }
            dcio::create_directory(beam_sidecar_path.str());
        }
    } else {
        WLOG(m_workerid, "already seen beam: " << beam_uid);
    }
}



void W2FWorker_Central::write_beamlet(BaseW2FWorker::W2FData& data) {
    std::unique_lock<std::mutex> lock(s_central_mutex);
    // create beamlet group
    std::ostringstream groupname;
    groupname << SparsifyManager::beams_data_group() + "/beam_" << std::setfill('0') << std::setw(5) << data.beam_header.beam_uid
              << std::setfill('\0')
              << "/" << "beamlet_" << std::setfill('0') << std::setw(5) << data.beamlet_header.beamlet_uid;

    WLOG(m_workerid, "creating group: " << "\"" << groupname.str() << "\"");
    auto beamletgroup = parent.m_h5file.createGroup(groupname.str());
    WLOG(m_workerid, "done creating group");
    WLOG(m_workerid, "Writing beamlet metadata");
    _write_beamlet_metadata(beamletgroup, data.beamlet_header, data.sparsedata);
    WLOG(m_workerid, "writing beamlet to file");
    _write_beamlet_to_hdf5(beamletgroup, data.sparsedata, data.beamlet_header);
    WLOG(m_workerid, "done writing beamlet to file");
}

// void W2FWorker_PerBeam::write_beamlet(BaseW2FWorker::W2FData& data) {
    /* write each beamlet to separate per-beam file named after beam_uid all in the subfolder
     * */
// }

void W2FWorker_PerBeamlet::write_beamlet(BaseW2FWorker::W2FData& data) {
    /* write each beamlet to separate file named after beam_uid and beamlet_uid all in the subfolder
     * */
    std::ostringstream targetfilename;
    targetfilename << "beam" << std::setfill('0') << std::setw(5) << data.beam_header.beam_uid <<
        std::setfill('\0') << "_beamlet" << std::setfill('0') << std::setw(5) << data.beamlet_header.beamlet_uid << std::setfill('\0') << ".bmlt";

    std::ostringstream sourcegroupname;
    sourcegroupname << SparsifyManager::beams_data_group() + "/beam_" << std::setfill('0') << std::setw(5) << data.beam_header.beam_uid
        << std::setfill('\0')
        << "/" << "beamlet_" << std::setfill('0') << std::setw(5) << data.beamlet_header.beamlet_uid;

    // open beamlet sidecar file
    std::ostringstream subpath;
    subpath << "beam_" << std::setfill('0') << std::setw(5) << data.beam_header.beam_uid << "/" + targetfilename.str();
    std::string filepath = parent.m_sidecar_path + subpath.str();
    {
        WLOG(m_workerid, "opening beamlet sidecar file: "<<filepath);
        auto h5file = H5::H5File(filepath, H5F_ACC_TRUNC);
        // targetgroup = targetgroup.createGroup(targetgroupname);
        auto targetgroup = h5file.openGroup("/");
        _write_beamlet_metadata(targetgroup, data.beamlet_header, data.sparsedata);
        _write_beamlet_to_hdf5(targetgroup, data.sparsedata, data.beamlet_header);
    }
    WLOG(m_workerid, "done writing beamlet to file");

    // create external link to beamlet file
    std::unique_lock<std::mutex> lock(s_central_mutex);
    WLOG(m_workerid, "creating group: " << "\"" << sourcegroupname.str() << "\"");
    auto sourcegroup = parent.m_h5file.createGroup(sourcegroupname.str());
    WLOG(m_workerid, "done creating group");

    WLOG(m_workerid, "Writing beamlet metadata");
    _write_beamlet_metadata(sourcegroup, data.beamlet_header, data.sparsedata);

    {
        std::string targetgroupname = "/lindex";
        std::string linkname = "lindex";
        WLOG(m_workerid, "creating link in "<<sourcegroupname.str()<<linkname<< " to: "<<filepath<<":"<<targetgroupname);
        H5Lcreate_external((dcio::get_basename(parent.m_sidecar_path) + "/" + subpath.str()).c_str(), targetgroupname.c_str(), sourcegroup.getId(), linkname.c_str(), H5P_DEFAULT, H5P_DEFAULT);
        WLOG(m_workerid, "done creating link");
    } {
        std::string targetgroupname = "/coeffs";
        std::string linkname = "coeffs";
        WLOG(m_workerid, "creating link in "<<sourcegroupname.str()<<linkname<< " to: "<<filepath<<":"<<targetgroupname);
        H5Lcreate_external((dcio::get_basename(parent.m_sidecar_path) + "/" + subpath.str()).c_str(), targetgroupname.c_str(), sourcegroup.getId(), linkname.c_str(), H5P_DEFAULT, H5P_DEFAULT);
        WLOG(m_workerid, "done creating link");
    }
}

void SparsifyManager::init(const std::string& h5filename, HEADER_PATIENT& patient_header, unsigned int num_srworkers, unsigned int num_w2fworkers, float sparsity_thresh, WRITE_STRATEGY write_strategy) {
    if (m_activated) { throw std::logic_error("Cannot initialize after activation of SparsifyManager\n"); }

    m_thresh = sparsity_thresh;
    m_write_strategy = write_strategy;

    if (write_strategy == WRITE_STRATEGY::CENTRAL) {
        num_w2fworkers = 1;
    }

    // open central hdf5 file to use for life of the worker
    try {
        m_h5file = H5::H5File(h5filename, H5F_ACC_TRUNC);
    } catch (H5::FileIException e) {
        std::cout << "Failed to open beamlet dose output file for writing. Please confirm that \""<<h5filename<<"\" is in a writable directory" << std::endl;
        //std::cout << e.getDetailMsg() << std::endl;
        throw;
    }

    if (write_strategy==WRITE_STRATEGY::PER_BEAMLET)
    {
        m_sidecar_path = h5filename + "_external/";
        LOG("Creating directory: " << m_sidecar_path);
        if (dcio::dir_exists(m_sidecar_path)) {
            dcio::remove_directory(m_sidecar_path, false, true);
        }
        dcio::create_directory(m_sidecar_path);
    }

    // write patient metadata to root group once
    H5::Group rootgroup = m_h5file.openGroup("/");
    _write_file_version(rootgroup, FTMAGIC, FTVERSIONMAJOR, FTVERSIONMINOR);
    _write_patient_metadata(rootgroup, patient_header);

    // create link hierarchy
    rootgroup.createGroup(beams_group());
    rootgroup.createGroup(beams_meta_group());
    rootgroup.createGroup(beams_data_group());

    //TODO: if nthreads==1 - use inline processing?
    // instantiate workers
    for (unsigned int ii=0; ii<num_srworkers; ++ii) {
        m_srworkers.emplace_back(new SRWorker(*this));
    }
    for (unsigned int ii=0; ii<num_w2fworkers; ++ii) {
        if (write_strategy == WRITE_STRATEGY::CENTRAL) {
            m_w2fworkers.emplace_back(new W2FWorker_Central(*this));
        } else if (write_strategy == WRITE_STRATEGY::PER_BEAMLET) {
            m_w2fworkers.emplace_back(new W2FWorker_PerBeamlet(*this));
        }
    }
    m_initialized = true;
}
void SparsifyManager::activate() {
    if (!m_initialized) { throw std::logic_error("Cannot activate SparsifyManager that hasn't been initialized\n"); }

    for (auto& srworker : m_srworkers) {
        srworker->start();  /* Launch thread */
    }
    for (auto& w2fworker : m_w2fworkers) {
        w2fworker->start();  /* Launch thread */
    }
}
void SparsifyManager::deactivate() {
    /* deactivating the manger will ensure that all items in each queue are processed (wqueue::deactivate())
     *   note: deactivate is a blocking call that only returns when wqueue is empty and "poison" signal
     *   has been sent to all workers. Thus it is safe to kill worker threads after deactivate returns
     * Once all items are processed for a queue, all attached workers will automatically die, then
     *   they will be joinable() */
    m_srqueue.deactivate();
    for (auto& srworker : m_srworkers) { srworker->join(); }
    m_w2fqueue.deactivate();
    for (auto& w2fworker : m_w2fworkers) { w2fworker->join(); }
}


unsigned int SRWorker::s_workerid = 0;
unsigned int BaseW2FWorker::s_workerid = 0;
std::mutex BaseW2FWorker::s_central_mutex{};
std::set<unsigned short> BaseW2FWorker::s_beams_seen{};
