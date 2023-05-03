#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

namespace wk
{
    class TrackingAction : public G4UserTrackingAction
    {
    public:
        TrackingAction();
        virtual ~TrackingAction(){};
    };
}

#endif