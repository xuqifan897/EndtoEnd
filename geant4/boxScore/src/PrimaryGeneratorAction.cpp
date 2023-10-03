#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorAction.h"
#include "argparse.h"
#include "PhantomDef.h"

bs::PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    this->fParticleGun = new G4ParticleGun(1);

    float energy = (*vm)["Energy"].as<float>() * MeV;
    this->SAD = (*vm)["SAD"].as<float>() * cm;
    this->beamletSize = (*vm)["beamlet-size"].as<float>() * cm;
    float SizeZ = 0.;
    for (int i=0; i<bs::GD->layers.size(); i++)
        SizeZ += std::get<1>(bs::GD->layers[i]);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    std::string particleName;
    G4ParticleDefinition* particle 
        = particleTable->FindParticle(particleName="gamma");
    this->fParticleGun->SetParticleDefinition(particle);
    this->fParticleGun->SetParticleEnergy(energy);
    this->fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., SizeZ - this->SAD));
}

bs::PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete this->fParticleGun;
}

void bs::PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // momentum sampling
    float isoplaneX = this->beamletSize * (G4UniformRand() - 0.5) * 2;
    float isoplaneY = this->beamletSize * (G4UniformRand() - 0.5) * 2;
    this->fParticleGun->SetParticlePosition(G4ThreeVector(isoplaneX, isoplaneY, this->SAD));
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}