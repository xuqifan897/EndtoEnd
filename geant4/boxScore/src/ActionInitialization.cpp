#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"

void bs::ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction(this->LocalScore));
}

void bs::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction(this->LocalScore));
}