#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

namespace bs
{
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        PrimaryGeneratorAction();
        virtual ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event*) override;
    
    private:
        G4ParticleGun* fParticleGun;
        float beamletSize;  // half value
        float SAD;
    };
}

#endif