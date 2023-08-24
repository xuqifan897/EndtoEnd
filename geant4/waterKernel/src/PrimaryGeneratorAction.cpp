#include "PrimaryGeneratorAction.h"
#include "argparse.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

wk::PrimaryGeneratorAction::PrimaryGeneratorAction()
    :fParticleGun(0)
{
    G4int n_particle = 1;
    G4ParticleGun* particleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName("gamma");
    G4ParticleDefinition* particle
        = particleTable->FindParticle(particleName);
    particleGun->SetParticleDefinition(particle);
    // particle momentum along the +z direction
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    
    G4float energy = getArg<float>("energy");
    particleGun->SetParticleEnergy(energy * MeV);
    
    G4float posZ = getArg<float>("posZ");
    particleGun->SetParticlePosition(G4ThreeVector(0., 0., posZ*cm));

    this->fParticleGun = particleGun;
}

#include "G4AutoLock.hh"
namespace wk
{
  G4Mutex PrimGenDestrMutex = G4MUTEX_INITIALIZER;
  G4Mutex PrimGenMutex = G4MUTEX_INITIALIZER;
}

wk::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    G4AutoLock lock(&wk::PrimGenDestrMutex);
    delete this->fParticleGun;
}

void wk::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}