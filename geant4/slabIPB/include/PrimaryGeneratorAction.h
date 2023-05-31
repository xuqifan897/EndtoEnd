#ifndef siPrimaryGeneratorAction_h
#define siPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4VPrimaryGenerator;
class G4Event;

namespace si
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