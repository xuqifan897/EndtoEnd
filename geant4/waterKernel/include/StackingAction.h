#ifndef StackingAction_h
#define StackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;
class G4VHitsCollection;

namespace wk
{
    class StackingAction : public G4UserStackingAction
    {
    public:
        StackingAction();
        virtual ~StackingAction();

        virtual G4ClassificationOfNewTrack 
            ClassifyNewTrack(const G4Track* aTrack);
    
    private:
        G4bool debugTrackingStacking;
    };
}

#endif