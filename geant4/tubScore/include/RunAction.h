#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

namespace ts
{
    class RunAction : public G4UserRunAction
    {
    public:
        RunAction() = default;
        ~RunAction() override = default;

        virtual G4Run* GenerateRun();
        void BeginOfRunAction(const G4Run*) override;
        void EndOfRunAction(const G4Run*) override;
    };
}

#endif