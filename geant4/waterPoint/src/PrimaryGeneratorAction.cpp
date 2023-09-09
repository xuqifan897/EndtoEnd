#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

#include "PrimaryGeneratorAction.h"
#include "argparse.h"
#include "config.h"
#include "EventInfo.h"

wp::PrimaryGeneratorAction::PrimaryGeneratorAction()
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

    // calculate the source position.
    int sourceOffsetUnitless = (*vm)["PhantomSZ"].as<int>() - (*vm)["PhantomDimZ"].as<int>() / 2;
    float sourceOffset = sourceOffsetUnitless * (*vm)["resolution"].as<float>() * cm * 2;
    this->fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., sourceOffset));
}

wp::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete this->fParticleGun;
}

void wp::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
    anEvent->SetUserInformation(new EventInfo());
}