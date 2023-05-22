#include "StackingAction.h"
#include "argparse.h"
#include "Trajectory.h"

wk::StackingAction::StackingAction()
    :G4UserStackingAction()
{
    this->debugTrackingStacking = getArg<bool>("debugTrackingStacking");
}

wk::StackingAction::~StackingAction()
{}

G4ClassificationOfNewTrack wk::StackingAction
    ::ClassifyNewTrack(const G4Track* aTrack)
{
    if (this->debugTrackingStacking)
    {
        G4cout << "ClassifyNewTrack processing Track " << 
            aTrack->GetTrackID() << G4endl;
    }
    G4ClassificationOfNewTrack classification = fUrgent;
    return classification;
}