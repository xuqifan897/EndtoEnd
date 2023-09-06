#include "Trajectory.h"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"

G4ThreadLocal G4Allocator<wp::Trajectory> * wp::myTrajectoryAllocator = 0;

wp::Trajectory::Trajectory(const G4Track* aTrack)
:G4VTrajectory()
{
    this->fParticleDefinition = aTrack->GetDefinition();
    this->fParticleName = this->fParticleDefinition->GetParticleName();
    this->fPDGCharge = this->fParticleDefinition->GetPDGCharge();
    this->fPDGEncoding = this->fParticleDefinition->GetPDGEncoding();
    if(fParticleName=="unknown")
    {
        G4PrimaryParticle*pp = aTrack->GetDynamicParticle()->GetPrimaryParticle();
        if(pp)
        {
            if(pp->GetCharge()<DBL_MAX) fPDGCharge = pp->GetCharge();
            this->fPDGEncoding = pp->GetPDGcode();
            if(pp->GetG4code()!=0)
            {
                this->fParticleName += " : ";
                this->fParticleName += pp->GetG4code()->GetParticleName();
            }
        }
    }
    this->fTrackID = aTrack->GetTrackID();
    this->fParentID = aTrack->GetParentID();
    this->fPositionRecord = new TrajectoryPointContainer();
    this->fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
    this->fMomentum = aTrack->GetMomentum();
    this->fVertexPosition = aTrack->GetPosition();
    this->fGlobalTime = aTrack->GetGlobalTime();
}

wp::Trajectory::~Trajectory()
{
    for (int i=0; i<this->fPositionRecord->size(); i++)
        delete (*(this->fPositionRecord))[i];
    this->fPositionRecord->clear();
    delete this->fPositionRecord;
}

void wp::Trajectory::ShowTrajectory(std::ostream& os) const
{
    os << G4endl << "TrackID = " << this->fTrackID
        << ": ParentID = " << this->fParentID << G4endl;
    os << "Particle name : " << this->fParticleName << " PDG code : " << this->fPDGEncoding
        << " Charge : " << this->fPDGCharge << G4endl;
    os << "Original momentum : " << G4BestUnit(this->fMomentum, "Energy") << G4endl;
    os << "Vertex : " << G4BestUnit(this->fVertexPosition, "Length")
        << " Global time: " << G4BestUnit(this->fGlobalTime, "Time") << G4endl;
    os << " Current trajectory has " << this->fPositionRecord->size() << " points." << G4endl;

    for (int i=0; i<this->fPositionRecord->size(); i++)
    {
        G4TrajectoryPoint* aTrajectoryPoint = 
            (G4TrajectoryPoint*)((*this->fPositionRecord)[i]);
        os << "Point[" << i << "]" 
            << " Position = " << aTrajectoryPoint->GetPosition() << G4endl;
    }
}

const std::map<G4String,G4AttDef>* wp::Trajectory::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01Trajectory",isNew);
    if (isNew)
    {
        G4String id("ID");
        (*store)[id] = G4AttDef(id,"Track ID","Bookkeeping","","G4int");

        G4String pid("PID");
        (*store)[pid] = G4AttDef(pid,"Parent ID","Bookkeeping","","G4int");

        G4String pn("PN");
        (*store)[pn] = G4AttDef(pn,"Particle Name","Bookkeeping","","G4String");

        G4String ch("Ch");
        (*store)[ch] = G4AttDef(ch,"Charge","Physics","e+","G4double");

        G4String pdg("PDG");
        (*store)[pdg] = G4AttDef(pdg,"PDG Encoding","Bookkeeping","","G4int");

        G4String imom("IMom");
        (*store)[imom] = G4AttDef(imom, "Momentum of track at start of trajectory",
                                    "Physics","G4BestUnit","G4ThreeVector");

        G4String imag("IMag");
        (*store)[imag] = 
            G4AttDef(imag, "Magnitude of momentum of track at start of trajectory",
                    "Physics","G4BestUnit","G4double");

        G4String vtxPos("VtxPos");
        (*store)[vtxPos] = G4AttDef(vtxPos, "Vertex position",
                                    "Physics","G4BestUnit","G4ThreeVector");

        G4String ntp("NTP");
        (*store)[ntp] = G4AttDef(ntp,"No. of points","Bookkeeping","","G4int");
    }
    return store;
}

std::vector<G4AttValue>* wp::Trajectory::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

    values->push_back
        (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

    values->push_back
        (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

    values->push_back(G4AttValue("PN",fParticleName,""));

    values->push_back
        (G4AttValue("Ch",G4UIcommand::ConvertToString(fPDGCharge),""));

    values->push_back
        (G4AttValue("PDG",G4UIcommand::ConvertToString(fPDGEncoding),""));

    values->push_back
        (G4AttValue("IMom",G4BestUnit(fMomentum,"Energy"),""));

    values->push_back
        (G4AttValue("IMag",G4BestUnit(fMomentum.mag(),"Energy"),""));

    values->push_back
        (G4AttValue("VtxPos",G4BestUnit(fVertexPosition,"Length"),""));

    values->push_back
        (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

    return values;
}

void wp::Trajectory::AppendStep(const G4Step* aStep)
{
    this->fPositionRecord->push_back(
        new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()));
}

void wp::Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
    if (! secondTrajectory) return;

    Trajectory* seco = (Trajectory*)secondTrajectory;
    for (int i=1; i<seco->GetPointEntries(); i++)
    {
        this->fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
    }
    delete (*seco->fPositionRecord)[0];
    seco->fPositionRecord->clear();
}