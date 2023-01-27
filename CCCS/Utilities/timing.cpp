#include "timing.h"
#include <cstdio>

/* PRIVATE */
double Timer::_difftime_ms(timepoint_t start, timepoint_t stop) {
    if (stop < start) {
        timepoint_t temp = start;
        start = stop;
        stop = temp;
    }
    return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000);
}


/* PUBLIC */
void Timer::reset() {
    m_start = hires_clock_t::now();
    counter = 0.0;
}
void Timer::start() {
    m_start = hires_clock_t::now();
    m_running = true;
}
void Timer::stop() {
    if (m_running) {
        m_stop = hires_clock_t::now();
        m_running = false;
        counter += _difftime_ms(m_start, m_stop);
    }
}
double Timer::time_elapsed() {
    // return time elapsed since start() in ms
    if (m_running) {
        // not stopped yet, keep running and print current time elapsed
        return counter + _difftime_ms(m_start, hires_clock_t::now());
    } else {
        return counter;
    }
}
void Timer::print_time_elapsed(const char* memo) {
    printf("** Time elapsed during \"%s\": ", memo);

    double _time_elapsed = time_elapsed();
    if (true || _time_elapsed < PRINT_MS) {
        printf("%.3f msec **\n", memo, _time_elapsed);
    } else if (_time_elapsed < PRINT_SEC) {
        printf("%.3f sec **\n", memo, _time_elapsed/1000);
    } else if (_time_elapsed < PRINT_MIN) {
        printf("%.3f min **\n", memo, _time_elapsed/1000/60);
    } else {
        printf("%.3f hrs **\n", memo, _time_elapsed/1000/60/60);
    }
}
void Timer::print_time_elapsed(const std::string& memo) {
    print_time_elapsed(memo.c_str());
}
void Timer::stop_print_time_elapsed(const char* memo) {
    stop();
    print_time_elapsed(memo);
}
void Timer::stop_print_time_elapsed(const std::string& memo) {
    stop_print_time_elapsed(memo.c_str());
}
void Timer::reset_print_time_elapsed(const char* memo) {
    print_time_elapsed(memo);
    reset();
}
void Timer::reset_print_time_elapsed(const std::string& memo) {
    reset_print_time_elapsed(memo.c_str());
}
void Timer::restart_print_time_elapsed(const char* memo) {
    print_time_elapsed(memo);
    start();
}
void Timer::restart_print_time_elapsed(const std::string& memo) {
    restart_print_time_elapsed(memo.c_str());
}
