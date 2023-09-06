#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include <atomic>

namespace wp
{
    extern std::atomic<bool> found;

    class EventAction : public G4UserEventAction
    {
    public:
        EventAction() = default;
        ~EventAction() = default;

        virtual void EndOfEventAction(const G4Event*);
    };
}

#endif