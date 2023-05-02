#include "TrackerSD.h"
#include "TrackerHit.h"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHandle.hh"
#include "G4ios.hh"

wk::TrackerSD::TrackerSD(G4String name)
    :G4VSensitiveDetector(name), fTrackerCollection(0)
{
    G4String HCname;
    collectionName.insert(HCname="trackerCollection");
}

wk::TrackerSD::~TrackerSD(){}

void wk::TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
    static int HCID = -1;
    this->fTrackerCollection = new TrackerHitsCollection(
        this->SensitiveDetectorName, collectionName[0]);
    if (HCID < 0)
        HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID, this->fTrackerCollection);
}

G4bool wk::TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    TrackerHit* newHit = new TrackerHit();
    newHit->SetEdep(edep);
    newHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
    newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
    fTrackerCollection->insert(newHit);
}

void wk::TrackerSD::EndOfEvent(G4HCofThisEvent*) {}