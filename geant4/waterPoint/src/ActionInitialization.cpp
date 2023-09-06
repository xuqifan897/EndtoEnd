#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "EventAction.h"
#include "TrackingAction.h"

void wp::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void wp::ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
    SetUserAction(new TrackingAction);
}