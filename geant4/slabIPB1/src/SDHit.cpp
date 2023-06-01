#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

#include "SDHit.h"

G4ThreadLocal G4Allocator<si::SDHit>* si::SDHitsAllocator = nullptr;

si::SDHit::SDHit()
    :G4VHit(), fEdep(0), fPos(0), fTrackID(-1) {}

si::SDHit::~SDHit() {}

void si::SDHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager)
    {
        G4Circle circle(fPos);
        circle.SetScreenSize(0.04);
        circle.SetFillStyle(G4Circle::filled);
        G4Color color(1.0, 0., 0.);
        G4VisAttributes attribs(color);
        circle.SetVisAttributes(attribs);
    }
}

void si::SDHit::Print()
{
    G4cout << "TrackID " << this->fTrackID << "   Position " << this->fPos 
        << "      :" << this->fEdep / keV << " [keV]" << G4endl;
}

const std::map<G4String, G4AttDef>* si::SDHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String, G4AttDef>* store 
        = G4AttDefStore::GetInstance("SDHit", isNew);
    if (isNew)
    {
        G4String hitType("HitType");
        (*store)[hitType] = G4AttDef(hitType, "Hit type", "Physics", "", "G4String");

        G4String trackID("trackID");
        (*store)[trackID] = G4AttDef(trackID, "Track ID", "Physics", "", "G4int");

        G4String energy("Energy");
        (*store)[energy] = G4AttDef(energy, "Energy Deposited", "Physics", 
            "G4BestUnit", "G4double");
        
        G4String pos("Pos");
        (*store)[pos] = G4AttDef(pos, "Position", "Physics", 
            "G4BestUnit", "G4ThreeVector");
    }
    return store;
}

std::vector<G4AttValue>* si::SDHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("HitType", "TrackerHit", ""));

    values->push_back(G4AttValue(
        "TrackID", G4UIcommand::ConvertToString(this->fTrackID), ""));
    
    values->push_back(G4AttValue("Energy", G4BestUnit(this->fEdep, "Energy"), ""));

    values->push_back(G4AttValue("Pos", G4BestUnit(this->fPos, "Length"), ""));

    return values;
}