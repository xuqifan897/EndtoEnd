#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

namespace wk
{
    class SteppingAction : public G4UserSteppingAction
    {
    public:
        SteppingAction();
        virtual ~SteppingAction();
    };
}

#endif