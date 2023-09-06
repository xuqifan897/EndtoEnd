#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"

#include "EventAction.h"
#include "argparse.h"
#include "Trajectory.h"

std::atomic<bool> wp::found(false);

void wp::EventAction::EndOfEventAction(const G4Event* evt)
{
    if ((*vm)["PrintTrajectory"].as<bool>())
    {
        G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
        int n_trajectories = 0;
        if (trajectoryContainer)
            n_trajectories = trajectoryContainer->entries();

        if (n_trajectories != 1)
        {
            bool expected = false;
            bool newValue = true;
            if (found.compare_exchange_strong(expected, newValue))
                for (int i=0; i<n_trajectories; i++)
                {
                    auto* trj = (*trajectoryContainer)[i];
                    trj->ShowTrajectory();
                }
        }
    }
}