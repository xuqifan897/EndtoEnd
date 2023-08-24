#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"


void ts::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void ts::ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
}