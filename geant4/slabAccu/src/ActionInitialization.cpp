#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"

void sa::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void sa::ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
}