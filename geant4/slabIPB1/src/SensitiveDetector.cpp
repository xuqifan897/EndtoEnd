#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHandle.hh"
#include "G4ios.hh"

#include "SensitiveDetector.h"
#include "SDHit.h"

si::SensitiveDetector::SensitiveDetector(G4String name, G4int index)
    :G4VSensitiveDetector(name), fSDHitsCollection(0), idx(index)
{
    // the parameter name is in the format "${layerID}${material}"
    G4String HCname(name + G4String("Collection"));
    collectionName.insert(HCname);
}

void si::SensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    this->fSDHitsCollection = new SDHitsCollection(
        this->SensitiveDetectorName, collectionName[0]);
    static int HCID = -1;
    if (HCID < 0)
        HCID = GetCollectionID(0);
    // for debug purposes
    G4cout << "HCID = " << HCID << G4endl;
    HCE->AddHitsCollection(HCID, this->fSDHitsCollection);
}

G4bool si::SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0) return false;

    SDHit* newHit = new SDHit();
    newHit->SetEdep(edep);
    newHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
    newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
    this->fSDHitsCollection->insert(newHit);
    return true;
}