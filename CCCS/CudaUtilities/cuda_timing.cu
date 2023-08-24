#include "cuda_timing.cuh"
#include <cstdio>

/* PRIVATE */
void CudaTimer::_init_events() {
    if (!m_init) {
        checkCudaErrors( cudaEventCreate(&m_start) );
        checkCudaErrors( cudaEventCreate(&m_stop) );
        m_init = true;
    }
}
float CudaTimer::_difftime_ms(cudaEvent_t start, cudaEvent_t stop) {
    if (stop < start) {
        cudaEvent_t temp = start;
        start = stop;
        stop = temp;
    }
    float elapsed_time_ms;
    checkCudaErrors( cudaEventSynchronize(m_stop) );
    checkCudaErrors( cudaEventElapsedTime(&elapsed_time_ms, m_start, m_stop) );
    return elapsed_time_ms;
}


/* PUBLIC */
void CudaTimer::reset() {
    m_started = m_stopped = false;
}
void CudaTimer::start() {
    _init_events();
    reset();
    checkCudaErrors( cudaEventRecord(m_start, m_stream) );
    m_started = true;
}
void CudaTimer::restart() { start(); }
void CudaTimer::stop() {
    if (m_started) {
        checkCudaErrors( cudaEventRecord(m_stop, m_stream) );
        m_stopped = true;
    }
}
float CudaTimer::time_elapsed() {
    // return time elapsed since start() in ms
    if (!m_started) { return 0.0f; }
    else if (!m_stopped) {
        // not stopped yet, keep running and print current time elapsed
        checkCudaErrors( cudaEventRecord(m_stop, m_stream) );
    }
    return _difftime_ms(m_start, m_stop);
}
void CudaTimer::print_time_elapsed(const char* memo) {
    float _time_elapsed = time_elapsed();
    printf("-- CUDA Time elapsed during \"%s\": ", memo);

    if (_time_elapsed < PRINT_MS) {
        printf("%.3f msec", memo, _time_elapsed);
    } else if (_time_elapsed < PRINT_SEC) {
        printf("%.3f sec", memo, _time_elapsed/1000);
    } else if (_time_elapsed < PRINT_MIN) {
        printf("%.3f min", memo, _time_elapsed/1000/60);
    } else {
        printf("%.3f hrs", memo, _time_elapsed/1000/60/60);
    }

    printf(" --\n");
}
void CudaTimer::print_time_elapsed(const std::string& memo) {
    print_time_elapsed(memo.c_str());
}
void CudaTimer::stop_print_time_elapsed(const char* memo) {
    stop();
    print_time_elapsed(memo);
}
void CudaTimer::stop_print_time_elapsed(const std::string& memo) {
    stop_print_time_elapsed(memo.c_str());
}
void CudaTimer::reset_print_time_elapsed(const char* memo) {
    print_time_elapsed(memo);
    reset();
}
void CudaTimer::reset_print_time_elapsed(const std::string& memo) {
    reset_print_time_elapsed(memo.c_str());
}
void CudaTimer::restart_print_time_elapsed(const char* memo) {
    print_time_elapsed(memo);
    start();
}
void CudaTimer::restart_print_time_elapsed(const std::string& memo) {
    restart_print_time_elapsed(memo.c_str());
}
