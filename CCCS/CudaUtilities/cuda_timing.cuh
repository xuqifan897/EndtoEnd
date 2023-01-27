#ifndef __CUDA_TIMING_H__
#define __CUDA_TIMING_H__

#include <helper_cuda.h>
#include <string>

// timing print format thresholds
#define PRINT_MS  1000
#define PRINT_SEC 1000 * 60
#define PRINT_MIN  1000 * 3600

class CudaTimer {
    private:
        bool m_init = false;
        cudaStream_t m_stream = 0;
        cudaEvent_t m_start = 0;
        cudaEvent_t m_stop = 0;
        bool m_started = false;
        bool m_stopped = false;

        void _init_events();
        float _difftime_ms(cudaEvent_t start, cudaEvent_t stop);

    public:
        CudaTimer() {}
        CudaTimer(cudaStream_t stream) : m_stream(stream) {}
        ~CudaTimer() {
            if (m_init) {
                checkCudaErrors( cudaEventDestroy(m_start) );
                checkCudaErrors( cudaEventDestroy(m_stop) );
            }
        }

        void reset();
        void start();
        void restart(); // alias for start()
        void stop();
        float time_elapsed();

        void print_time_elapsed(const char* memo);
        void print_time_elapsed(const std::string& memo);
        void stop_print_time_elapsed(const char* memo);
        void stop_print_time_elapsed(const std::string& memo);
        void reset_print_time_elapsed(const char* memo);
        void reset_print_time_elapsed(const std::string& memo);
        void restart_print_time_elapsed(const char* memo);
        void restart_print_time_elapsed(const std::string& memo);
};

class CudaAutoTimer : public CudaTimer {
    // derived version of Timer which automatically starts upon instantiation
    public:
        CudaAutoTimer() { start(); }
        CudaAutoTimer(cudaStream_t stream) : CudaTimer(stream) {
            start();
        }
};

#endif // __CUDA_TIMING_H__
