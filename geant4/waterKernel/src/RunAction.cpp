#include "RunAction.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

wk::RunAction::RunAction()
    :G4UserRunAction()
{}

wk::RunAction::~RunAction()
{}

void wk::RunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << "start." << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

void wk::RunAction::EndOfRunAction(const G4Run*)
{}