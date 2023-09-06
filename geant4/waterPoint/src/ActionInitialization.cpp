#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"

void wp::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void wp::ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
}