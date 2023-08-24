#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "PrimaryGeneratorAction.h"
#include "argparse.h"
#include "PhantomDef.h"

ts::PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    G4int n_particle = 1;
    this->fParticleGun = new G4ParticleGun(n_particle);

    // default particle kinematic
    G4float energy = (*vm)["Energy"].as<float>() * MeV;
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    std::string particleName;
    G4ParticleDefinition* particle
        = particleTable->FindParticle(particleName="gamma");
    this->fParticleGun->SetParticleDefinition(particle);
    this->fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    this->fParticleGun->SetParticleEnergy(energy);
    this->fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -GD->sizeZ));
}

ts::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete this->fParticleGun;
}

void ts::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}