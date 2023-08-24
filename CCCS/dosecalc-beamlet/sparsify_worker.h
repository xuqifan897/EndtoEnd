#ifndef __SPARSIFY_WORKER_H__
#define __SPARSIFY_WORKER_H__

#include <string>
#include <unordered_set>

#include "Utilities/thread.h"
#include "Utilities/wqueue.h"
#include "DoseCalcIO/sparse_io.h"
#include "DoseCalcIO/macros.h"
#include "H5Cpp.h" // namespace H5::

// Data storage to be consumed by worker thread
struct SparsifyData {
    std::string    filename;
    const float*   dense_array = nullptr;
    ArrayProps     props;
    HEADER_BEAMLET beamlet_header;
    HEADER_BEAM    beam_header;
    ROIMaskList* roi_list = nullptr;
};

const extern float _default_sparsity_threshold;
const extern int _default_compress_lvl;
const extern int _default_chunksize;
class SparsifyWorkerBase : public Thread {
    protected:
        static unsigned int s_workerid;

        wqueue<SparsifyData>&    m_queue;
        unsigned int             m_workerid;
        std::unordered_set<uint> m_beams_seen;
        const float              m_thresh = _default_sparsity_threshold;
        int                m_compress_lvl = _default_compress_lvl;
        int                m_chunksize = _default_chunksize;

        void _assign_workerid() { m_workerid = s_workerid++;}

    public:
        SparsifyWorkerBase(wqueue<SparsifyData>& queue, const float sparsity_thresh=_default_sparsity_threshold)
            : m_queue{queue}, m_thresh{sparsity_thresh}
        { _assign_workerid(); }
        virtual ~SparsifyWorkerBase();

        void set_compress_lvl(int c) { m_compress_lvl = c; }
        int get_compress_lvl() { return m_compress_lvl; }
        void set_chunksize(int c) { m_chunksize = c; }
        int get_chunksize() { return m_chunksize; }

        uint workerid() { return m_workerid; }
        virtual void join();
        virtual void run() {};
};
class BeamletSparsifyWorker : public SparsifyWorkerBase {
    public:
        BeamletSparsifyWorker(wqueue<SparsifyData>& queue, const float sparsity_thresh=_default_sparsity_threshold)
            : SparsifyWorkerBase(queue, sparsity_thresh) {}

        void run();
};

class FullSparsifyWorker : public SparsifyWorkerBase {
    protected:
        std::string m_h5filename;
        H5::H5File m_h5file;

        int write_sparse_full(H5::H5File&, SparsifyData& sdata);

    public:
        HEADER_PATIENT* patient_header;

        FullSparsifyWorker(
                wqueue<SparsifyData>& queue,
                const std::string h5filename,
                HEADER_PATIENT* patient_header,
                const float sparsity_thresh=_default_sparsity_threshold
        );
        void run();
};

#endif // __SPARSIFYWORKER_H__
