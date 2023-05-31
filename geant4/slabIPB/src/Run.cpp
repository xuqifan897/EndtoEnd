#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4Threading.hh"

#include "Run.h"
#include "PhantomDef.h"
#include "argparse.h"
#include "SDHit.h"

si::Run::Run():
    G4Run(),
    sizeX(GD->sizeX), sizeY(GD->sizeY), sizeZ(GD->sizeZ),
    fullSizeX(2 * GD->sizeX), fullSizeY(2 * GD->sizeY), fullSizeZ(2 * GD->sizeZ),
    resX(GD->resX), resY(GD->resY), resZ(GD->resZ),
    fullResX(2 * GD->resX), fullResY(2 * GD->resY), fullResZ(2 * GD->resZ),
    dimX(GD->dimX), dimY(GD->dimY), dimZ(GD->dimZ)
{
    if (GD == nullptr)
    {
        G4cout << "Material geometry is not initialized, "
            "please initialize it by calling \"si::GD = new "
            "si::GeomDef();\" before the detector construction" << G4endl;
        exit(1);
    }
    this->recordEventLog = getArg<bool>("recordEventLog");
    G4SDManager * SDman = G4SDManager::GetSDMpointer();

    this->HitsCollectionIDs = std::vector<G4int>();
    this->Names = std::vector<G4String>();
    G4int layerID = 1;
    for (auto it = GD->layers.begin(); it != GD->layers.end(); it++)
    {
        G4String& material = std::get<0>(*it);
        G4String collectionName = std::to_string(layerID) 
            + material + G4String("Collection");
        layerID ++;
        this->HitsCollectionIDs.push_back(
            SDman->GetCollectionID(collectionName));
        this->Names.push_back(collectionName);
    }

    // initialize fHitsMap. The scoring region is 
    // assumed to be the same as the phanotm
    this->fHitsMap = std::vector<G4double>(this->dimX * this->dimY * this->dimZ);

    if (G4Threading::IsMasterThread())
    {
        for (G4int i=0; i<GD->layers.size(); i++)
            G4cout << "Layer " << this->Names[i] << ", collection ID: " 
                << this->HitsCollectionIDs[i] << G4endl;
    }
}

si::Run::~Run() {}

void si::Run::RecordEvent(const G4Event* evt)
{
    this->numberOfEvent ++;
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    if (! HCE)
        return;

    SDHitsCollection * THC;
    // iterate over all detectors
    for (int i=0; i<this->HitsCollectionIDs.size(); i++)
    {
        THC = 0;
        THC = (SDHitsCollection*)(HCE->GetHC(this->HitsCollectionIDs[i]));
        if (!THC)
            continue;
        int n_hit = THC->entries();
        if (this->recordEventLog)
        {
            G4cout << "Event ID: " << evt->GetEventID() << ". Collection ID:" 
            << i+1 << ". Number of hits: " << n_hit << G4endl;
        }
        if (n_hit == 0)
            continue;
        for (int j=0; j<n_hit; j++)
        {
            SDHit* aHit = (*THC)[j];
            const G4ThreeVector& aHitPos = aHit->GetPos();
            size_t idx = this->index(aHitPos);
            this->fHitsMap[idx] += aHit->GetEdep();
            if (this->recordEventLog)
            {
                G4cout << "(TrackID, HitID): (" << aHit->getTrackID() << ", " 
                    << j << "). Hit Position: (" << aHitPos << "). Index: " 
                    << idx << ". Detector: " << this->Names[i] << G4endl;;
            }
        }
    }
}

size_t si::Run::index(G4ThreeVector aHitPos)
{
    aHitPos[0] += this->sizeX;
    aHitPos[1] += this->sizeY;
    aHitPos[2] += this->sizeZ;
    size_t idX = int(aHitPos[0] / this->fullResX);
    size_t idY = int(aHitPos[1] / this->fullResY);
    size_t idZ = int(aHitPos[2] / this->fullResZ);
    size_t result = ((idX * this->dimY) + idY) * this->dimZ + idZ;
    return result;
}

void si::Run::Merge(const G4Run* aRun)
{
    const si::Run* localRun = static_cast<const si::Run*>(aRun);
    this->numberOfEvent += localRun->numberOfEvent;
    for (int i=0; i<this->fHitsMap.size(); i++)
        this->fHitsMap[i] += localRun->fHitsMap[i];
}

void si::Run::writeHitsMap(G4String path) const
{
    std::ofstream file(path);
    if (file.is_open())
    {
        file.write((char*)(this->fHitsMap.data()),
            this->fHitsMap.size()*sizeof(G4double));
        file.close();
        G4cout << "Array data written successfully" << G4endl;
    }
    else
        G4cout << "Unable to open file: " << path << G4endl;
}