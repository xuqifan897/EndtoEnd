#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4VPrimaryGenerator;
class G4Event;

namespace wk
{
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        PrimaryGeneratorAction();
        virtual ~PrimaryGeneratorAction();
    
        virtual void GeneratePrimaries(G4Event* anEvent);
    
    private:
        G4VPrimaryGenerator* fParticleGun;
    };
}

#endif