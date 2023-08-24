#ifndef __WQUEUE_H__
#define __WQUEUE_H__

#include <iostream>
#include <stdexcept>
#include <queue>
#include <mutex>
#include <condition_variable>


#if 0
    #define LOG(message)     std::cout << message << std::endl
    #define WLOG(wid,message) LOG("w" << wid << ":  " << message)
#else
    #define LOG(message) {}
    #define WLOG(wid,message) {}
#endif


// defines a thread-compatible FIFO work queue that can be fed by producer thread and consumed by worker threads
// Items are enqueued using wqueue::push
// Items are "popped" using try_pop, which returns bool indicating success of pop and whether
// valid data now occupies the parameter reference "item" or if it is left untouched
//
// When the queue will no longer be added to, but items remain and further processing must occur,
//   the caller should execute wqueue::deactivate() which will change the state of the object,
//   prompting it to finish processing remaining items, then send kill signals to attached
//   workers when all items are removed. This is done by having workers monitor the return bool
//   of try_pop for a false signal then breaking the try_pop endless loop and exiting the
//   ::run() method the manager of the workers is then responsible for "join()"-ing all managed
//   worker threads which will block until all threads are killed
template<typename T>
class wqueue {
    private:
        std::queue<T> m_queue;
        std::recursive_mutex m_mutex;
        std::condition_variable_any m_condvar;
        bool m_activated;
        bool m_killnow;

    public:
        wqueue() : m_activated{true}, m_killnow(false) {}

        // push by move
        void push(T&& item) {
            std::unique_lock<std::recursive_mutex> lock(m_mutex);
            if (m_activated) {
                LOG("Pushing onto the stack");
                m_queue.push(std::move(item));
                LOG("current wqueue size: " << size());
                LOG("done pushing, notifying one");
                m_condvar.notify_one();
                LOG("done notifying");
            } else {
                throw std::logic_error("wqueue::push() was attempted on a deactivated queue");
            }
        }
        // push by copy
        void push(T& item) {
            std::unique_lock<std::recursive_mutex> lock(m_mutex);
            if (m_activated) {
                LOG("Pushing onto the stack");
                m_queue.push(item);
                LOG("current wqueue size: " << size());
                LOG("done pushing, notifying one");
                m_condvar.notify_one();
                LOG("done notifying");
            } else {
                throw std::logic_error("wqueue::push() was attempted on a deactivated queue");
            }
        }

        bool try_pop(T& item, int workerid=-1) {
            WLOG(workerid, "Entered try_pop()");
            std::unique_lock<std::recursive_mutex> lock(m_mutex);

            WLOG(workerid, "current wqueue size: " << size());
            while (size() <= 0 && !m_killnow) {
                WLOG(workerid, "waiting to dequeue");
                m_condvar.wait(lock);
                /* -- Wait for signal -- */
                WLOG(workerid, "done waiting");
            }

            if (size() > 0) {
                WLOG(workerid, "dequeuing");
                item = std::move(m_queue.front());
                m_queue.pop();
                WLOG(workerid, "done dequeuing - " << &item);
                return true;
            } else { // causes calling thread to terminate/join
                WLOG(workerid, "Kill switch triggered in try_pop()");
                return false;
            }
        }

        int size() {
            std::unique_lock<std::recursive_mutex> lock(m_mutex);
            return m_queue.size();
        }

        void deactivate() {
            // DO NOT USE LOCK HERE -- BLOCKS AND CAUSES ENDLESS LOOP
            m_activated = false;
            while (size() > 0) {
                // finish remaining tasks
                LOG("firing finish-up notification");
                LOG("current wqueue size: " << size());
                m_condvar.notify_one();
            }

            // kill all consumer threads (signaled with size==0)
            LOG("firing kill notification to all");
            m_killnow = true;
            m_condvar.notify_all();
        }
        void reset() {
            std::unique_lock<std::recursive_mutex> lock(m_mutex);
            LOG("current wqueue size: " << size());
            LOG("clearing queue");
            m_queue = std::queue<T>(); // clear queue
            LOG("current wqueue size: " << size());
            m_activated = true;
            m_killnow = false;
            m_condvar.notify_all();
        }

};

#endif // __WQUEUE_H__
