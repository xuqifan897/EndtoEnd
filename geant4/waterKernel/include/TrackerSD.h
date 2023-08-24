#ifndef wkTrackerSD_h
#define wkTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.h"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

namespace wk
{

    class TrackerSD : public G4VSensitiveDetector
    {
    public:
        TrackerSD(G4String name);
        virtual ~TrackerSD();

        virtual void Initialize(G4HCofThisEvent* HCE);
        virtual void EndOfEvent(G4HCofThisEvent* HCE);
    
    protected:
        virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    
    private:
        TrackerHitsCollection * fTrackerCollection;
    };

}

#endif