#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "argparse.h"
#include "PhantomDef.h"
#include "PrimaryGeneratorAction.h"

si::PrimaryGeneratorAction::PrimaryGeneratorAction()
{   
    G4int n_particle = 1;
    this->fParticleGun = new G4ParticleGun(n_particle);

    // default particle kinematic
    G4float energy = si::getArg<float>("Energy") * MeV;
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle
        = particleTable->FindParticle(particleName="gamma");
    this->fParticleGun->SetParticleDefinition(particle);
    this->fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    this->fParticleGun->SetParticleEnergy(energy);
    this->fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -GD->sizeZ));
}

si::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete this->fParticleGun;
}

void si::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}