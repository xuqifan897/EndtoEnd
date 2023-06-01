#ifndef siPrimaryGeneratorAction_h
#define siPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

namespace si
{
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        PrimaryGeneratorAction();
        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event*) override;
    
    private:
        G4ParticleGun* fParticleGun;
    };
}

#endif