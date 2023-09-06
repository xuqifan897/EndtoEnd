#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

namespace wp
{

    class TrackingAction : public G4UserTrackingAction
    {
    public:
        TrackingAction();
        virtual ~TrackingAction() = default;

        virtual void PreUserTrackingAction(const G4Track*);
    };
}

#endif