#include <string>
#include <iostream>
#include <atomic>

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "G4THitsMap.hh"
#include "G4Event.hh"

#include "Run.h"
#include "argparse.h"
#include "config.h"
#include "EventInfo.h"
#include "EventAction.h"

std::atomic<int> wp::globalCount = 0;

wp::Run::Run()
{
    // // firstly, initialize hit collection ID
    this->logFrequency = (*vm)["logFrequency"].as<int>();
    this->nParticles = (*vm)["nParticles"].as<int>();
    this->PhantomDimXY = (*vm)["PhantomDimXY"].as<int>();
    this->PhantomDimZ = (*vm)["PhantomDimZ"].as<int>();
    this->resolution = (*vm)["resolution"].as<float>() * cm;
    this->HalfPhantomSizeZ = this->PhantomDimZ * this->resolution;
    int numSenseDet = this->PhantomDimXY * this->PhantomDimZ;

    this->HitsMaps = std::vector<std::tuple<std::string, int, G4THitsMap<G4double>*>>(numSenseDet);
    G4SDManager* pSDman = G4SDManager::GetSDMpointer();
    for (int i=0; i<numSenseDet; i++)
    {
        std::string name = std::string("det") + std::to_string(i+1) + std::string("/Edep");
        int id = pSDman->GetCollectionID(name);
        std::get<0>(this->HitsMaps[i]) = name;
        std::get<1>(this->HitsMaps[i]) = id;
        std::get<2>(this->HitsMaps[i]) = nullptr;
    }

    // The order is the same, z, y, x
    this->PhantomSZ = (*vm)["PhantomSZ"].as<int>();
    this->PhantomBottom = (*vm)["PhantomBottom"].as<int>();
    this->marginZ = this->HalfPhantomSizeZ - this->PhantomBottom * 2 * this->resolution;
    int kernelDimZ = this->PhantomSZ + this->PhantomBottom;
    int kernelDimXY = (*vm)["PhantomDimXY"].as<int>();
    this->kernel = std::vector<std::vector<std::vector<double>>>(
        kernelDimZ, std::vector<std::vector<double>>(kernelDimXY, 
        std::vector<double>(kernelDimXY)));
}

wp::Run::~Run()
{}

void wp::Run::RecordEvent(const G4Event* anEvent)
{
    G4VUserEventInformation* VUserInfo = anEvent->GetUserInformation();
    EventInfo* UserInfo = static_cast<EventInfo*>(VUserInfo);
    bool IfInteraction = UserInfo->GetIfInteraction();
    double InteractionZ = UserInfo->GetInitInteractionZ();
    if ((! IfInteraction) || InteractionZ >= this->marginZ)
        return;  // If no interaction happened, just return
    this->numberOfEvent ++;

    int newCount = globalCount.fetch_add(1);
    if (newCount % this->logFrequency == 0)
        std::cout << "Event number: " << newCount 
            << " / " << this->nParticles << std::endl;

    // determine the range of the HitsMap that should be registered
    int InteractionDim = int((InteractionZ + this->HalfPhantomSizeZ) 
        / (2 * this->resolution));
    int sliceBegin = InteractionDim - this->PhantomSZ;
    int sliceEnd = InteractionDim + this->PhantomBottom;

    for (int i=0; i<this->PhantomSZ+this->PhantomBottom; i++)
    {
        int idxZ = sliceBegin + i;
        int detIDZ = idxZ * this->PhantomDimXY;
        for (int idxY=0; idxY<this->PhantomDimXY; idxY++)
        {
            int detID = detIDZ + idxY;
            int HCID = std::get<1>(this->HitsMaps[detID]);
            auto hitsCollection = static_cast<G4THitsMap<G4double>*>(
                anEvent->GetHCofThisEvent()->GetHC(HCID));
            const auto & map = *(hitsCollection->GetMap());
            for (const auto & it : map)
                this->kernel[i][idxY][it.first] += *(it.second);
        }
    }
}

void wp::Run::Merge(const G4Run* aRun)
{
    const Run* bRun = static_cast<const Run*>(aRun);

    // Order: z, y, x
    for (int i=0; i<this->PhantomSZ + this->PhantomBottom; i++)
        for (int j=0; j<this->PhantomDimXY; j++)
            for (int k=0; k<this->PhantomDimXY; k++)
                this->kernel[i][j][k] += bRun->kernel[i][j][k];
    this->numberOfEvent += bRun->GetNumberOfEvent();
}