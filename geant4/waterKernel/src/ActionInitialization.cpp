#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "EventAction.h"
#include "StackingAction.h"
#include "TrackingAction.h"

wk::ActionInitialization::ActionInitialization(){}

wk::ActionInitialization::~ActionInitialization(){}

void wk::ActionInitialization::Build() const
{
    //
    SetUserAction(new wk::PrimaryGeneratorAction);
    //
    SetUserAction(new wk::RunAction);
    //
    SetUserAction(new wk::EventAction);

    //
    SetUserAction(new wk::StackingAction);
    //
    SetUserAction(new wk::TrackingAction);

}

void wk::ActionInitialization::BuildForMaster() const
{
    G4UserRunAction* run_action = new RunAction;
    SetUserAction(run_action);
}