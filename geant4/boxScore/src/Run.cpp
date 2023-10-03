#include <string>
#include <iostream>

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"

#include "Run.h"
#include "PhantomDef.h"
#include "argparse.h"

bs::Run::Run()
{
    int dimXY = (*bs::vm)["dimXY"].as<int>();
    int SegZ = (*bs::vm)["SegZ"].as<int>();
    this->HitsMaps.reserve(dimXY*SegZ);
    for (int i=0; i<dimXY*SegZ; i++)
    {
        std::string name = "SD" + std::to_string(i+1) + "/Edep";
        int id = G4SDManager::GetSDMpointer()->GetCollectionID(name);
        this->HitsMaps.push_back(std::make_tuple(name, id, new G4THitsMap<G4double>()));
    }
}

bs::Run::~Run()
{
    for (int i=0; i<this->HitsMaps.size(); i++)
        delete std::get<2>(this->HitsMaps[i]);
}

void bs::Run::RecordEvent(const G4Event* anEvent)
{
    for (int i=0; i<this->HitsMaps.size(); i++)
    {
        int HCID = std::get<1>(this->HitsMaps[i]);
        auto hitsCollection = static_cast<G4THitsMap<G4double>*>(
            anEvent->GetHCofThisEvent()->GetHC(HCID));
        *std::get<2>(this->HitsMaps[i]) += *hitsCollection;
    }
}

void bs::Run::Merge(const G4Run* aRun)
{
    const Run* bRun = static_cast<const Run*>(aRun);
    for (int i=0; i<this->HitsMaps.size(); i++)
        *std::get<2>(this->HitsMaps[i]) += *std::get<2>(bRun->HitsMaps[i]);
}