#ifndef wkTrackerHit_h
#define wkTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include "config.h"

class G4AttDef;
class G4AttValue;

namespace wk
{

#if PHASE == 0
    class TrackerHit : public G4VHit
    {
    public:
        TrackerHit();
        virtual ~TrackerHit();

        inline void *operator new(size_t);
        inline void operator delete(void *aHit);

        virtual void Draw();
        virtual void Print();
        virtual const std::map<G4String, G4AttDef>* GetAttDefs() const;
        virtual std::vector<G4AttValue>* CreateAttValues() const;
    
    public:
        inline void SetEdep(G4double de)
        { this->fEdep = de; }
        inline void AddEdep(G4double de)
        { this->fEdep += de; }
        inline G4double GetEdep() const
        { return this->fEdep; }
        inline void SetPos(G4ThreeVector xyz)
        { this->fPos = xyz;}
        inline void SetTrackID(G4int i)
        { this->fTrackID = i; }
        inline G4int GetTrackID() const
        { return this->fTrackID; }
    
    private:
        G4double fEdep;
        G4ThreeVector fPos;
        G4int fTrackID;

    };
#endif

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

// I don't know what the following object is, just copy the example RE01
extern G4ThreadLocal G4Allocator<TrackerHit> * TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
    if(!TrackerHitAllocator)
        TrackerHitAllocator = new G4Allocator<TrackerHit>;
    return (void *) TrackerHitAllocator->MallocSingle();
}

inline void TrackerHit::operator delete(void *aHit)
{
    TrackerHitAllocator->FreeSingle((TrackerHit*) aHit);
}

}

#endif