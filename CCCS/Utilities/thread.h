#ifndef __THREAD_WRAPPER_H__
#define __THREAD_WRAPPER_H__

#include <thread>

// defines base class for any threaded execution
// derived must define run() method which will be run once after thread is launched
// Thread is launched by executing Thread::start()
// Thread is killed by executing Thread::join()
class Thread {
    protected:
        std::thread self;
        bool m_running, m_detached;

    public:
        Thread() : m_running{false}, m_detached{false} {}
        virtual ~Thread() {}

        std::thread::id get_id();
        void start();
        void join();
        void detach();

        // abstract method - must redefine in derived class
        virtual void run() = 0;
};

#endif // __THREAD_WRAPPER_H__
