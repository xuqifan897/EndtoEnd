#include "TrackingAction.h"
#include "Trajectory.h"
#include "argparse.h"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

wk::TrackingAction::TrackingAction()
    :G4UserTrackingAction()
{
    this->debugTrackingStacking = getArg<bool>("debugTrackingStacking");
}

void wk::TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    this->fpTrackingManager->SetStoreTrajectory(true);
    this->fpTrackingManager->SetTrajectory(new Trajectory(aTrack));
}

void wk::TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries)
    {
        size_t nSeco = secondaries->size();
        if (nSeco > 0 && this->debugTrackingStacking)
        {
            for (size_t i=0; i<nSeco; i++)
            {
                G4cout << "PostUserTrackingAction processing track: " 
                    << ((*secondaries)[i])->GetTrackID() << G4endl;
            }
        }
    }
}