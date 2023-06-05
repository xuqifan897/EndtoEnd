#ifndef siSDHit_h
#define siSDHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4AttDef;
class G4AttValue;

namespace si
{
    class SDHit : public G4VHit
    {
    public:
        SDHit();
        virtual ~SDHit();

        inline void * operator new(size_t);
        inline void operator delete(void* aHit);

        virtual void Draw();
        virtual void Print();
        virtual const std::map<G4String, G4AttDef>* GetAttDefs() const;
        virtual std::vector<G4AttValue>* CreateAttValues() const;
    
    public:
        inline void SetEdep(G4double de)
            {this->fEdep = de;}
        inline G4double& GetEdep()
            {return this->fEdep;}
        inline void SetPos(G4ThreeVector xyz)
            {this->fPos = xyz;}
        inline G4ThreeVector& GetPos()
            {return this->fPos;}
        inline void SetTrackID( G4int i)
            {this->fTrackID = i;}
        inline G4int& getTrackID()
            {return this->fTrackID;}
    private:
        G4double fEdep;
        G4ThreeVector fPos;
        G4int fTrackID;
    };

    typedef G4THitsCollection<SDHit> SDHitsCollection;

    // I don't understand what the following object is, just copy the example RE01
    extern G4ThreadLocal G4Allocator<SDHit>* SDHitsAllocator;

    inline void* SDHit::operator new(size_t)
    {
        if(!SDHitsAllocator)
            SDHitsAllocator = new G4Allocator<SDHit>;
        return (void*) SDHitsAllocator->MallocSingle();
    }

    inline void SDHit::operator delete(void* aHit)
    {
        SDHitsAllocator->FreeSingle((SDHit*) aHit);
    }
}

#endif