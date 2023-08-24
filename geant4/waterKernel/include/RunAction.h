#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "config.h"

class G4Run;

namespace wk
{
    class RunAction : public G4UserRunAction
    {
    public:
        RunAction();
        virtual ~RunAction();

    #if PHASE > 0
        virtual G4Run* GenerateRun();
    #endif
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);
    };
}

#endif