#include <string>
#include <iostream>
#include <atomic>

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"

#include "Run.h"
#include "argparse.h"

std::atomic<int> wp::globalCount = 0;

wp::Run::Run()
{
    // // firstly, initialize hit collection ID
    this->logFrequency = (*vm)["logFrequency"].as<int>();
    this->nParticles = (*vm)["nParticles"].as<int>();
    int PhantomDimXY = (*vm)["PhantomDimXY"].as<int>();
    int PhantomDimZ = (*vm)["PhantomDimZ"].as<int>();
    int numSenseDet = PhantomDimXY * PhantomDimZ;
    this->HitsMaps = std::vector<std::tuple<std::string, int, G4THitsMap<G4double>*>>(numSenseDet);
    G4SDManager* pSDman = G4SDManager::GetSDMpointer();
    for (int i=0; i<numSenseDet; i++)
    {
        std::string name = std::string("det") + std::to_string(i+1) + std::string("/Edep");
        int id = pSDman->GetCollectionID(name);
        std::get<0>(this->HitsMaps[i]) = name;
        std::get<1>(this->HitsMaps[i]) = id;
        std::get<2>(this->HitsMaps[i]) = new G4THitsMap<double>();
    }
}

wp::Run::~Run()
{
    // delete std::get<2>(this->HitsMap);

    for (int i=0; i<this->HitsMaps.size(); i++)
    {
        delete std::get<2>(this->HitsMaps[i]);
    }
}

void wp::Run::RecordEvent(const G4Event* anEvent)
{
    for (int i=0; i<this->HitsMaps.size(); i++)
    {
        auto hitsCollection = static_cast<G4THitsMap<G4double>*>(
            anEvent->GetHCofThisEvent()->GetHC(std::get<1>(this->HitsMaps[i])));
        *(std::get<2>(this->HitsMaps[i])) += *hitsCollection;
    }
    this->numberOfEvent ++;

    // for logging purposes
    int newCount = globalCount.fetch_add(1);
    if (newCount % this->logFrequency == 0)
        std::cout << "Event number: " << newCount 
            << " / " << this->nParticles << std::endl;
}

void wp::Run::Merge(const G4Run* aRun)
{
        const Run* bRun = static_cast<const Run*>(aRun);
    for (int i=0; i<this->HitsMaps.size(); i++)
    {
        *std::get<2>(this->HitsMaps[i]) += *std::get<2>(bRun->HitsMaps[i]);
    }
    this->numberOfEvent += bRun->GetNumberOfEvent();
}