#include "EventAction.h"
#include "TrackerHit.h"
#include "Trajectory.h"
#include "config.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

wk::EventAction::EventAction()
    :G4UserEventAction(), fTrackerCollID(-1)
{}

wk::EventAction::~EventAction()
{}

void wk::EventAction::BeginOfEventAction(const G4Event*)
{
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    if (this->fTrackerCollID<0)
    {
        G4String colNam;
        this->fTrackerCollID = SDman->GetCollectionID(colNam="trackerCollection");
    }
}

void wk::EventAction::EndOfEventAction(const G4Event* evt)
{
    G4cout << ">>> Summary of Event " << evt->GetEventID() << G4endl;
    if(evt->GetNumberOfPrimaryVertex()==0)
    {
        G4cout << "Event is empty." << G4endl;
        return;
    }

    if (this->fTrackerCollID < 0) return;

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    TrackerHitsCollection* THC = 0;
    if (HCE)
    {
        THC = (TrackerHitsCollection*)(HCE->GetHC(this->fTrackerCollID));
    }

    #if PHASE == 0  // print every hit
    if (THC)
    {
        int n_hit = THC->entries();
        G4cout << G4endl;
        G4cout << "Tracker hits " << 
            "--------------------------------------------------------------"
            << G4endl;
        G4cout << n_hit << " hits are stored in TrackerHitsCollection."
            << G4endl;
        G4cout << "List of hits in tracker" << G4endl;
        for (int i=0; i<n_hit; i++)
        {
            (*THC)[i]->Print();
        }
    }
    #endif

    #if PHASE == 0 // print primary particles
    G4cout << G4endl;
    G4cout << "Primary particles " <<
        "--------------------------------------------------------------"
        << G4endl;
    G4int n_vertex = evt->GetNumberOfPrimaryVertex();
    for (G4int iv=0; iv<n_vertex;iv++)
    {
        G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
        G4cout << G4endl;
        G4cout << "Primary vertex "
            << G4ThreeVector(pv->GetX0(), pv->GetY0(), pv->GetZ0())
            << "   at t = " << (pv->GetT0())/ns << " [ns]" << G4endl;
        
        if (true)
        {
            G4PrimaryParticle* pp = pv->GetPrimary();
            while(pp)
            {
                PrintPrimary(pp, 0);
                pp->GetNext();
            }
        }
    }
    #endif

    #if PHASE == 0
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    // extract the trajectories and print them
    G4cout << G4endl;
    G4cout << "Trajectories in tracker " <<
        "--------------------------------------------------------------"
        << G4endl;
    if (true)
    {
        for(G4int i=0; i<n_trajectories; i++)
        {
            Trajectory* trj = 
                (Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
            trj->ShowTrajectory();
        }
    }
    #endif
}

void wk::EventAction::PrintPrimary(G4PrimaryParticle* pp, G4int ind)
{
    for (G4int ii=0; ii<=ind; ii++)
        G4cout << "   ";
    G4cout << "==PDGcode " << pp->GetPDGcode() << " ";
    if (pp->GetG4code()!=0)
        G4cout << "(" << pp->GetG4code()->GetParticleName() << ")";
    else
        { G4cout << "is not defined in G4"; }
    G4cout << " " << pp->GetMomentum()/GeV << " [GeV] ";
    if(pp->GetTrackID()<0)
        { G4cout << G4endl; }
    else
        { G4cout << ">>> G4Track ID " << pp->GetTrackID() << G4endl; }

    G4PrimaryParticle* daughter = pp->GetDaughter();
    while(daughter)
    {
        PrintPrimary(daughter,ind+1);
        daughter = daughter->GetNext();
    }
}