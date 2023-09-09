#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>

#include "EventAction.h"
#include "argparse.h"
#include "Trajectory.h"
#include "config.h"
#include "EventInfo.h"
#include "TrajectoryPoint.h"

std::atomic<bool> wp::found(false);

float wp::EventAction::marginZ(0.);
float wp::EventAction::sourceOffset(0.);
wp::EventAction::EventAction()
{
    int PhantomDimZ = (*vm)["PhantomDimZ"].as<int>();
    int PhantomBottom = (*vm)["PhantomBottom"].as<int>();
    float resolution = (*vm)["resolution"].as<float>() * cm;
    float offset = 2 * resolution * ((float)PhantomDimZ / 2. - PhantomBottom);
    marginZ = offset;
    std::cout << "The limit in the z direction is " << 
        marginZ / cm << " cm." << std::endl;
    
    int sourceOffsetUnitless = (*vm)["PhantomSZ"].as<int>() - 
        (float)PhantomDimZ / 2;
    sourceOffset = sourceOffsetUnitless * 2 * resolution;
}

void wp::EventAction::EndOfEventAction(const G4Event* evt)
{
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    int n_trajectories = 0;
    if (trajectoryContainer)
        n_trajectories = trajectoryContainer->entries();

    if (n_trajectories == 1)
        return;

    // At first, I thought the first vertex point of a track should be 
    // the interaction point. But then I found that Some vertex point 
    // can be before the primary particle source. So I'll try another 
    // method to do so.
#if false
    for (int i=0; i<trajectoryContainer->size(); i++)
    {
        auto* VTrj = (*trajectoryContainer)[i];
        Trajectory* trj = static_cast<Trajectory*>(VTrj);
        if (trj->GetParentID()==1)
        {
            // This particle is generated by the primary interaction
            const G4ThreeVector & vertex = trj->GetVertexPosition();
            double interZ = vertex[2];

            // for debug purposes
            if (interZ < sourceOffset)
            {
                // Strange things happens. The initial interaction 
                // point is before the primary particle generator
                std::cerr << "The Z coordinate of the interaction is: " << interZ / cm 
                    << "cm, while the source Z coordinate is: " << sourceOffset / cm 
                    << "cm." << std::endl;
            }

            if (interZ < marginZ)
            {
                G4VUserEventInformation* VUserInfo = evt->GetUserInformation();
                EventInfo* UserInfo = static_cast<EventInfo*>(VUserInfo);
                UserInfo->SetIfInteraction(true);
                UserInfo->SetInitInteractionZ(interZ);
            }
            break;
        }
    }
#endif

    // Assume the first trajectory corresponds to the primary particle
    auto* VTrj = (*trajectoryContainer)[0];
    Trajectory* trj = static_cast<Trajectory*>(VTrj);

    // for debug purposes
    #if false
    if (trj->GetTrackID() != 1)
    {
        std::cerr << "The first trajectory does not "
            "correspond to the primary particle" << std::endl;
        exit(EXIT_FAILURE);
    }
    #endif

    int interactionIdx = 0;
    bool interaction_ = false;

    //  Here, we get the interaction point by analyzing the trajectory
#if false
    for (int i=0; i<trj->GetPointEntries(); i++)
    {
        const auto VPoint = trj->GetPoint(i);
        TrajectoryPoint* point = static_cast<TrajectoryPoint*>(VPoint);
        const G4ThreeVector & position = point->GetPosition();
        // assume the first point after the interaction 
        // should be deviated from the centerline
        if (std::abs(position[0])>EPS || std::abs(position[1]>EPS))
        {
            interactionIdx = i - 1;
            interaction_ = true;
            break;
        }
    }
#endif

#if true
    //  Here, we get the interaction point by analyzing momentum
    if (trj->GetFlag())
    {
        interactionIdx = trj->GetInterIdx();
        interaction_ = true;
    }
#endif

    //  If no interaction happened, or all trajectory 
    //  points are on the original direction, probably in the last interaction, 
    //  all energy is released to secondary particles.
    if(! interaction_)
    {
        const auto VLastPoint = trj->GetPoint(trj->GetPointEntries()-1);
        TrajectoryPoint* point = static_cast<TrajectoryPoint*>(VLastPoint);
        const G4ThreeVector & position = point->GetPosition();
        if (position[2] < marginZ)
        {
            //  The last interaction point was before the margin
            interaction_ = true;
            interactionIdx = trj->GetPointEntries()-1;
        }
    }

    // for debug purposes
    #if false
    if (interactionIdx < 0)
    {
        std::cerr << "The first interaction point is " 
            << interactionIdx << std::endl;
        exit(EXIT_FAILURE);
    }
    #endif

    if (interaction_)
    {
        G4VTrajectoryPoint* VInterPoint = trj->GetPoint(interactionIdx);
        TrajectoryPoint* InterPoint = static_cast<TrajectoryPoint*>(VInterPoint);
        double interZ = InterPoint->GetPosition()[2];

        // for debug purposes
        #if false
        if (interZ < sourceOffset)
        {
            // Strange things happens. The initial interaction 
            // point is before the primary particle generator
            std::cerr << "The Z coordinate of the interaction is: " << interZ / cm 
                << "cm, the source Z coordinate is: " << sourceOffset / cm 
                << "cm. The index of the trajectory point of interaction is " 
                << interactionIdx << std::endl;
            trj->ShowTrajectory();
            exit(EXIT_FAILURE);
        }
        #endif

        G4VUserEventInformation* VUserInfo = evt->GetUserInformation();
        EventInfo* UserInfo = static_cast<EventInfo*>(VUserInfo);
        UserInfo->SetIfInteraction(true);
        UserInfo->SetInitInteractionZ(interZ);
    }
}