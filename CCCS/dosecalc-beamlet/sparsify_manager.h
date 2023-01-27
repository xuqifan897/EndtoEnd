#ifndef __SPARSIFY_MANAGER_H__
#define __SPARSIFY_MANAGER_H__

#include <string>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <memory>
#include <mutex>

#include "DoseCalcIO/sparse_io.h"
#include "Utilities/thread.h"
#include "Utilities/wqueue.h"


class SparsifyManager; // fwd decl.
class ROIMaskList; // forward decl.

enum class SPARSIFY_STRATEGY {
    INLINE,
    THREADED,
};
enum class WRITE_STRATEGY {
    CENTRAL,
    PER_BEAM,
    PER_BEAMLET,
};


class SRWorker : public Thread {
    friend class SparsifyManager;
    private:
        SparsifyManager& parent;
        static unsigned int s_workerid;
        unsigned int m_workerid;
        void _assign_workerid() { m_workerid = s_workerid++;}
    public:
        SRWorker(SparsifyManager& parent) : parent{parent} { _assign_workerid(); }     /* default const. */
        SRWorker(SRWorker& other) : parent{other.parent} { _assign_workerid(); };  /* Copy const. */
        SRWorker(SRWorker&& other) : parent{other.parent} { _assign_workerid(); }; /* Move const. */
        ~SRWorker() { --s_workerid; }

        uint workerid() { return m_workerid; }
        void run();

        struct SRData {
            SRData() {};
            // move constructor - no copy because of std::unique_ptr
            SRData(SRData&& other) {
                props = std::move(other.props);
                beamlet_header = std::move(other.beamlet_header);
                beam_header = std::move(other.beam_header);
                roi_list = std::move(other.roi_list);
                dense_array = std::unique_ptr<float[]>(std::move(other.dense_array));
            }
            // move assign op.
            SRData& operator=(SRData&& other) {
                props = std::move(other.props);
                beamlet_header = std::move(other.beamlet_header);
                beam_header = std::move(other.beam_header);
                roi_list = std::move(other.roi_list);
                dense_array = std::unique_ptr<float[]>(std::move(other.dense_array));
                return *this;
            }

            // Data Members
            std::unique_ptr<float[]> dense_array;
            ArrayProps             props;
            HEADER_BEAMLET         beamlet_header;
            HEADER_BEAM            beam_header;
            ROIMaskList*           roi_list = nullptr;
        };
};

class BaseW2FWorker : public Thread {
    public:
        BaseW2FWorker(SparsifyManager& parent) : parent{parent} { _assign_workerid(); }
        BaseW2FWorker(const BaseW2FWorker& other) : parent{other.parent} { _assign_workerid(); }; /* Copy const. */

        SparsifyManager& parent;
        struct W2FData {
            SparseData     sparsedata;
            HEADER_BEAMLET beamlet_header;
            HEADER_BEAM    beam_header;
        };

        void run();

    protected:
        static unsigned int s_workerid;
        unsigned int m_workerid;
        void _assign_workerid() { m_workerid = s_workerid++;}

        static std::mutex s_central_mutex;
        static std::set<unsigned short> s_beams_seen;
        void write_beam_metadata(W2FData& data);
        virtual void write_beamlet(W2FData& data) = 0;
};

class W2FWorker_Central : public BaseW2FWorker {
    protected:
        void write_beamlet(BaseW2FWorker::W2FData& data);
    public:
        W2FWorker_Central(SparsifyManager& parent) : BaseW2FWorker(parent) {}
        W2FWorker_Central(const BaseW2FWorker& other) : BaseW2FWorker(other.parent) {}
};
// class W2FWorker_PerBeam : public BaseW2FWorker {
//     protected:
//         void write_beamlet(BaseW2FWorker::W2FData& data);

//     public:
//         W2FWorker_PerBeam(SparsifyManager& parent) : BaseW2FWorker(parent) {}
//         W2FWorker_PerBeam(const BaseW2FWorker& other) : BaseW2FWorker(other.parent) {}
// };
class W2FWorker_PerBeamlet : public BaseW2FWorker {
    protected:
        void write_beamlet(BaseW2FWorker::W2FData& data);

    public:
        W2FWorker_PerBeamlet(SparsifyManager& parent) : BaseW2FWorker(parent) {}
        W2FWorker_PerBeamlet(const BaseW2FWorker& other) : BaseW2FWorker(other.parent) {}
};


// manages a pool of workers and exposes wqueue-like interaction methods
const extern float _default_sparsity_threshold;
const extern int _default_compress_lvl;
const extern int _default_chunksize;
const WRITE_STRATEGY _default_write_strategy = WRITE_STRATEGY::CENTRAL;
class SparsifyManager {
    friend class BaseW2FWorker;
    friend class W2FWorker_Central;
    friend class W2FWorker_PerBeamlet;
    friend class SRWorker;
    private:
        H5::H5File m_h5file;
        std::string m_sidecar_path;
        bool m_activated=false;
        bool m_initialized=false;
        float m_thresh = _default_sparsity_threshold;
        WRITE_STRATEGY m_write_strategy;
        wqueue<SRWorker::SRData> m_srqueue {};
        wqueue<BaseW2FWorker::W2FData> m_w2fqueue {};
        std::vector<std::unique_ptr<SRWorker> > m_srworkers;
        std::vector<std::unique_ptr<BaseW2FWorker> > m_w2fworkers;

        // define standard paths in h5 file
        static std::string beams_group() { return "/beams"; }
        static std::string beams_meta_group() { return beams_group() + "/metadata"; }
        static std::string beams_data_group() { return beams_group() + "/data"; }

    public:
        std::mutex m_memblock_mutex;
        std::condition_variable m_memblock_condvar;
        std::stack<std::unique_ptr<float[]> > m_memblock_stack {};
        void init(const std::string& h5filename, HEADER_PATIENT& patient_header, unsigned int num_srworkers=1, unsigned int num_w2fworkers=1,
                float sparsity_thresh=_default_sparsity_threshold, WRITE_STRATEGY write_strategy=_default_write_strategy);
        SparsifyManager() {}
        SparsifyManager(const std::string& h5filename, HEADER_PATIENT& patient_header, unsigned int num_srworkers=1, unsigned int num_w2fworkers=1,
                float sparsity_thresh=_default_sparsity_threshold, WRITE_STRATEGY write_strategy=_default_write_strategy)
        { init(h5filename, patient_header, num_srworkers, num_w2fworkers, sparsity_thresh, write_strategy); }

        WRITE_STRATEGY get_write_strategy() { return m_write_strategy; }
        int get_num_srworkers() { return m_srworkers.size(); }
        int get_num_w2fworkers() { return m_w2fworkers.size(); }

        std::unique_ptr<float[]> get_memblock(bool blocking=true) {
            /* block until a memory address is available to return */
            LOG("get_memblock: enter - wait for non-empty stack" << std::endl);
            std::unique_lock<std::mutex> lock(m_memblock_mutex);
            if (!blocking && m_memblock_stack.empty()) {
                return nullptr;
            }
            while (m_memblock_stack.empty()) { m_memblock_condvar.wait(lock); }

            LOG("get_memblock: locking mutex (owns: "<<lock.owns_lock()<<")" << std::endl);
            std::unique_ptr<float[]> mem_addr = std::move(m_memblock_stack.top());
            m_memblock_stack.pop();
            LOG("get_memblock: got memblock, "<<m_memblock_stack.size()<<"left in stack unlocking mutex"<<std::endl);
            return mem_addr;
        }
        void push_memblock(std::unique_ptr<float[]>& addr) {
            LOG("push_memblock: enter - wait for mutex" << std::endl);
            std::unique_lock<std::mutex> lock(m_memblock_mutex);
            LOG("push_memblock: locking mutex (owns: "<<lock.owns_lock()<<")" << std::endl);
            m_memblock_stack.push(std::move(addr));
            m_memblock_condvar.notify_one();
            LOG("push_memblock: added memblock to stack. " << m_memblock_stack.size() << " available" << std::endl);
            LOG("push_memblock: exit; mutex unlock" << std::endl);
        }

        void activate();
        void push(SRWorker::SRData& item) { m_srqueue.push(std::move(item)); }
        int srqueue_size () { return m_srqueue.size(); }
        int w2fqueue_size () { return m_w2fqueue.size(); }
        void deactivate();
};


#endif //__SPARSIFY_MANAGER_H__
