#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"

void si::ActionInitialization::Build() const
{   
    //
    SetUserAction(new si::PrimaryGeneratorAction);
    //
    SetUserAction(new si::RunAction);
}

void si::ActionInitialization::BuildForMaster() const
{
    SetUserAction(new si::RunAction);
}