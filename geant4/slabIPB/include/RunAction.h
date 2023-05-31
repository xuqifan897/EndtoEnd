#ifndef siRunAction_h
#define siRunAction_h 1

#include "G4UserRunAction.hh"

class G4Run;

namespace si
{
    class RunAction : public G4UserRunAction
    {
    public:
        RunAction();
        virtual ~RunAction();

        virtual G4Run* GenerateRun();
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);
    };
}

#endif