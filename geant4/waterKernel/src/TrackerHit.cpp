#include "TrackerHit.h"

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

G4ThreadLocal G4Allocator<wk::TrackerHit> * wk::TrackerHitAllocator = 0;

wk::TrackerHit::TrackerHit()
    :G4VHit(), fEdep(0), fPos(0), fTrackID(-1) {}

wk::TrackerHit::~TrackerHit() {}

void wk::TrackerHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fPos);
        circle.SetScreenSize(0.04);
        circle.SetFillStyle(G4Circle::filled);
        G4Color colour(1., 0., 0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
    }
}

const std::map<G4String, G4AttDef>* wk::TrackerHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String, G4AttDef>* store 
        = G4AttDefStore::GetInstance("TrackerHit", isNew);
    if (isNew)
    {
        G4String hitType("HitType");
        (*store)[hitType] = G4AttDef(hitType, "Hit Type", "Physics", "", "G4String");

        G4String trackID("TrackID");
        (*store)[trackID] = G4AttDef(trackID, "Track ID", "Physics", "", "G4int");

        G4String energy("Energy");
        (*store)[energy] = G4AttDef(energy, "Energy Deposited", "Physics", 
            "G4BestUnit", "G4double");
        
        G4String pos("Pos");
        (*store)[pos] = G4AttDef(pos, "Position", 
            "Physics", "G4BestUnit", "G4ThreeVector");
    }
    return store;
}

std::vector<G4AttValue>* wk::TrackerHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("HitType", "TrackerHit", ""));

    values->push_back(
        G4AttValue("TrackID", G4UIcommand::ConvertToString(this->fTrackID), ""));
    
    values->push_back(G4AttValue("Energy", G4BestUnit(this->fEdep, "Energy"), ""));

    values->push_back(G4AttValue("Pos", G4BestUnit(this->fPos, "Length"), ""));

    return values;
}

void wk::TrackerHit::Print()
{
    G4cout << "TrackID " << this->fTrackID << "   Position " << this->fPos << "       : "
        << this->fEdep / keV << " [keV]" << G4endl;
}