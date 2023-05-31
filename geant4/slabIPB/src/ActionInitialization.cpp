#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"

si::ActionInitialization::ActionInitialization() {}

si::ActionInitialization::~ActionInitialization() {}

void si::ActionInitialization::Build() const
{
    //
    SetUserAction(new si::PrimaryGeneratorAction);
    //
    SetUserAction(new si::RunAction);
}

void si::ActionInitialization::BuildForMaster() const
{
    G4UserRunAction* run_action = new si::RunAction;
    SetUserAction(run_action);
}