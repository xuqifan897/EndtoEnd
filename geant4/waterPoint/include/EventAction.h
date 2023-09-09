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
        EventAction();
        ~EventAction() = default;

        virtual void EndOfEventAction(const G4Event*);
    
    private:
        static float marginZ;
        static float sourceOffset;
    };
}

#endif