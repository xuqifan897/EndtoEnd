#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "PrimaryGeneratorAction.h"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "argparse.h"
#include "PhantomDef.h"

si::PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    G4int n_particle = 1;
    G4ParticleGun* ParticleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName("gamma");
    G4ParticleDefinition* particle 
        = particleTable->FindParticle(particleName);
    ParticleGun->SetParticleDefinition(particle);
    
    // particle momentum along +z direction
    ParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    G4float energy = getArg<float>("Energy") * MeV;
    ParticleGun->SetParticleEnergy(energy);

    // the origin of particle is at the bottom of phantom in z direction
    ParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -GD->sizeZ));

    this->fParticleGun = ParticleGun;
}

#include "G4AutoLock.hh"

// I don't know what it means
namespace si
{
    G4Mutex PrimaryGenDestrMutex = G4MUTEX_INITIALIZER;
    G4Mutex PrimaryGenMutex = G4MUTEX_INITIALIZER;
}

si::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    G4AutoLock lock(&PrimaryGenDestrMutex);
    delete this->fParticleGun;
}

void si::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}