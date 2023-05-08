#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "config.h"

class G4PrimaryParticle;

namespace wk
{
    class EventAction : public G4UserEventAction
    {
    public:
        EventAction();
        virtual ~EventAction();

        virtual void BeginOfEventAction(const G4Event*);
    #if PHASE == 0
        virtual void EndOfEventAction(const G4Event*);
    #endif

    private:
        G4int fTrackerCollID;
        void PrintPrimary(G4PrimaryParticle* pp, G4int ind);
    };
}

#endif