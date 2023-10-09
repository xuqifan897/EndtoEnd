#include "DummyActionInitialization.h"
#include "PrimaryGeneratorAction.h"

void bs::DummyActionInitialization::BuildForMaster() const
{}

void bs::DummyActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
}