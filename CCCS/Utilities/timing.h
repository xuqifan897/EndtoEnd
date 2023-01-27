#ifndef __TIMING_H__
#define __TIMING_H__

#include <chrono>
#include <string>

// timing print format thresholds
#define PRINT_MS  1000
#define PRINT_SEC 1000 * 60
#define PRINT_MIN  1000 * 3600

typedef std::chrono::high_resolution_clock             hires_clock_t;
typedef std::chrono::high_resolution_clock::time_point timepoint_t;

class Timer {
    private:
        timepoint_t m_start;
        timepoint_t m_stop;
        bool m_running = false;

        double _difftime_ms(timepoint_t start, timepoint_t stop);
        double counter = 0.0;

    public:
        void reset();
        void start();
        void stop();
        double time_elapsed();

        void print_time_elapsed(const char* memo);
        void print_time_elapsed(const std::string& memo);
        void stop_print_time_elapsed(const char* memo);
        void stop_print_time_elapsed(const std::string& memo);
        void reset_print_time_elapsed(const char* memo);
        void reset_print_time_elapsed(const std::string& memo);
        void restart_print_time_elapsed(const char* memo);
        void restart_print_time_elapsed(const std::string& memo);
};

class AutoTimer : public Timer {
    // derived version of Timer which automatically starts upon instantiation
    public:
        AutoTimer() { start(); }
};

#endif // __TIMING_H__
