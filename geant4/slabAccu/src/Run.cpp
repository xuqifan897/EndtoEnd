#include <string>
#include <iostream>

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"

#include "Run.h"
#include "PhantomDef.h"

sa::Run::Run()
{
    // firstly, initialize hit collection ID
    for (int i=0; i<GD->layers.size(); i++)
    {
        std::string name = std::string("SD") + std::to_string(i+1) + 
            std::string("/Edep");
        int id = G4SDManager::GetSDMpointer()->GetCollectionID(name);
        this->HitsMaps.push_back(std::make_tuple(name, id, G4THitsMap<G4double>()));
    }

    // log
    if (G4Threading::IsMasterThread())
    {
        for (auto it=this->HitsMaps.begin(); it!=this->HitsMaps.end(); it++)
            std::cout << std::get<0>(*it) << ":" << std::setw(10) << 
                std::right << std::get<1>(*it) << std::endl;
    }
}

void sa::Run::RecordEvent(const G4Event* anEvent)
{
    for (auto it=this->HitsMaps.begin(); it!=this->HitsMaps.end(); it++)
    {
        int HCID = std::get<1>(*it);
        auto hitsCollection
            = static_cast<G4THitsMap<G4double>*>(
                anEvent->GetHCofThisEvent()->GetHC(HCID));
        std::get<2>(*it) += *hitsCollection;
    }
}

void sa::Run::Merge(const G4Run* aRun)
{
    Run* bRun = (Run*)aRun;
    for (int i=0; i<this->HitsMaps.size(); i++)
    {
        std::get<2>(this->HitsMaps[i]) += std::get<2>(bRun->HitsMaps[i]);
    }
}