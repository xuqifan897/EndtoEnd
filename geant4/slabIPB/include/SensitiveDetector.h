#ifndef siSensitiveDetector_h
#define siSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "SDHit.h"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

namespace si
{
    class SensitiveDetector : public G4VSensitiveDetector
    {
    public:
        SensitiveDetector(G4String name);
        SensitiveDetector(SensitiveDetector& old) = delete;
        SensitiveDetector(SensitiveDetector&& old) = delete;
        SensitiveDetector& operator=(const SensitiveDetector old) = delete;

        virtual ~SensitiveDetector();

        virtual void Initialize(G4HCofThisEvent* HCE);
        virtual void EndOfEvent(G4HCofThisEvent* HCE);
    
    protected:
        virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    
    private:
        SDHitsCollection* fSDHitsCollection;
    };
}

#endif