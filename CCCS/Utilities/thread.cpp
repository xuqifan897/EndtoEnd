#include "thread.h"

std::thread::id Thread::get_id() {
    return self.get_id();
}
void Thread::start() {
    self = std::thread(&Thread::run, this); // exception if failure
    m_running = true;
}
void Thread::join() {
    if (m_running && self.joinable()) { self.join(); }
    m_running = false;
}
void Thread::detach() {
    if (m_running && !m_detached && self.joinable()) { self.detach(); }
    m_detached = true;
}
