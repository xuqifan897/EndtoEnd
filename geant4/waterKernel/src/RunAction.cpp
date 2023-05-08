#include "RunAction.h"
#include "wkRun.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

wk::RunAction::RunAction()
    :G4UserRunAction()
{}

wk::RunAction::~RunAction()
{}

#if PHASE > 0
G4Run* wk::RunAction::GenerateRun()
{ return new wk::Run; }
#endif

void wk::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

void wk::RunAction::EndOfRunAction(const G4Run* aRun)
{
#if PHASE > 0
    // if (IsMaster())
    // {}
#endif
}