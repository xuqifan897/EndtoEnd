#include "TrackingAction.h"
#include "Trajectory.h"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

wk::TrackingAction::TrackingAction()
    :G4UserTrackingAction()
{}

void wk::TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    this->fpTrackingManager->SetStoreTrajectory(true);
    this->fpTrackingManager->SetTrajectory(new Trajectory(aTrack));
}