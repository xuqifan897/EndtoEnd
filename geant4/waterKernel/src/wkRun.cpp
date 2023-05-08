#include "wkRun.h"
#include "argparse.h"
#include "TrackerHit.h"
#include "Trajectory.h"

#include "globals.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include <vector>

wk::Run::Run()
    :G4Run(), numberOfEvent(0)
{
    G4String kernelType = getArg<std::string>("kernelType");
    if (kernelType == G4String("point"))
        this->pointKernel = true;
    else if (kernelType == G4String("IPB"))
        this->pointKernel = false;
    else
    {
        std::cerr << "unkown kernelType: " << kernelType;
        std::cerr << ". Assume the kernelType to be \"point\"";
        this->pointKernel = true;
    }

    this->RecordEventLog = getArg<bool>("recordEventLog");

    // get fTrackerCollID
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;
    this->fTrackerCollID = SDman->GetCollectionID(colNam="trackerCollection");
    
    if (this->pointKernel)
    {
        this->fSizeX = getArg<float>("kernelSizeX") * cm;
        this->fSizeY = getArg<float>("kernelSizeY") * cm;
        this->fSizeZ = getArg<float>("kernelSizeZ") * cm;

        this->fPosZ = getArg<float>("kernelPosZ") * cm;
    }
    else
    {
        this->fSizeX = getArg<float>("sizeX") * cm;
        this->fSizeY = getArg<float>("sizeY") * cm;
        this->fSizeZ = getArg<float>("sizeZ") * cm;

        // this->fPosZ = - this->fSizeZ;
    }

    this->fResX = getArg<float>("kernelResX") * cm;
    this->fResY = getArg<float>("kernelResY") * cm;
    this->fResZ = getArg<float>("kernelResZ") * cm;

    this->fFullResX = 2 * this->fResX;
    this->fFullResY = 2 * this->fResY;
    this->fFullResZ = 2 * this->fResZ;

    this->fDimX = int(this->fSizeX / this->fResX);
    this->fDimY = int(this->fSizeY / this->fResY);
    this->fDimZ = int(this->fSizeZ / this->fResZ);

    // then, recalculate the size
    this->fSizeX = this->fDimX * this->fResX;
    this->fSizeY = this->fDimY * this->fResY;
    this->fSizeZ = this->fDimZ * this->fResZ;

    if (! this->pointKernel)
        this->fPosZ = - this->fSizeZ;

    this->fFullSizeX = this->fSizeX * 2;
    this->fFullSizeY = this->fSizeY * 2;
    this->fFullSizeZ = this->fSizeZ * 2;

    bool kernelDimOdd = getArg<bool>("kernelDimOdd");
    if (kernelDimOdd)
    {
        // we only care about fDimX and fDimY
        if (this->fDimX % 2 == 0)
            this->fDimX -= 1;
        if (this->fDimY % 2 == 0)
            this->fDimY -= 1;
    }
    else
    {
        if (this->fDimX % 2 == 1)
            this->fDimX -= 1;
        if (this->fDimY % 2 == 1)
            this->fDimY -= 1;
    }

    if (this->pointKernel)
    {
        G4float kernelFrontMargin = this->fSizeZ - this->fPosZ;
        G4float kernelBackMargin = this->fSizeZ + this->fPosZ;
        G4float phantomSizeZ = getArg<float>("sizeZ") * cm;
        this->fValidRangeZ0 = - (phantomSizeZ - kernelBackMargin);
        this->fValidRangeZ1 = phantomSizeZ - kernelFrontMargin;
    }

    this->fHitsMap = std::vector<G4double>(this->fDimX * this->fDimY * this->fDimZ);
    this->logGeom();
}


wk::Run::~Run()
{}


size_t wk::Run::index(G4ThreeVector position, bool& validFlag, 
    int& hitID, G4double& edep)
{
    // This function transforms the position to the index of HitsMap
    // Firstly, extract the interaction point w.r.t. the base (0,0,0)
    // of fHitsMap. The parameters int& hitID and bool& validFlag are
    // for logging only.

    float posX, posY, posZ;
    posX = position[0];
    posY = position[1];
    posZ = position[2];

    // firstly, calculate the position w.r.t. the origin of the kernel
    posX += this->fSizeX;
    posY += this->fSizeY;
    posZ = this->fPosZ + posZ + this->fSizeZ;

    validFlag = (0 <= posX) && (posX < this->fFullSizeX)
        && (0 <= posY) && (posY < this->fFullSizeY)
        && (0 <= posZ) && (posZ < this->fFullSizeZ);

    if (! this->RecordEventLog)
        if (! validFlag)
            return 0;

    size_t idX, idY, idZ, result;
    idX = int(posX / this->fFullResX);
    idY = int(posY / this->fFullResY);
    idZ = int(posZ / this->fFullResZ);

    result = ((idX * this->fDimY) + idY) * this->fDimZ + idZ;

    if (this->RecordEventLog)
    {
        G4cout << "hit ID: " << hitID << G4endl;
        G4cout << "Original coordinates: (" << position[0] / cm 
            << ", " << position[1] / cm << ", " << position[2] / cm << ") [cm]." << G4endl;
        G4cout << "Position relative to base: (" << posX / cm 
            << ", " << posY / cm << ", " << posZ / cm << ") [cm]. "
            << "validFlag: " << validFlag <<  G4endl;
        G4cout << "coordinate relative to base: ( " << idX << 
            ", " << idY << ", " << idZ << "). 1D coordinate: "
            << result << G4endl;
        G4cout << "energy: " << edep / keV << "[keV].\n" << G4endl;
    }

    return result;
}


void wk::Run::RecordEvent(const G4Event* evt)
{
    G4double firstInteractionPointZ;
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    // Since the primary particle is gamma, the 
    // trajectory Container has at lease one entry
    Trajectory* trj = (Trajectory*)((*(evt->GetTrajectoryContainer()))[0]);
    if (this->pointKernel)
    {
        // if it is the point kernel, we firstly 
        // extract the first interaction point
        firstInteractionPointZ = trj->GetFirstInteractionPointZ();
        // to test whether the firstInteractionPoint is within the valid range
        if ((firstInteractionPointZ < this->fValidRangeZ0) || 
            (firstInteractionPointZ >= this->fValidRangeZ1))
        return;
    }
    else
        firstInteractionPointZ = trj->GetEntrancePointZ();

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    TrackerHitsCollection* THC = 0;
    if (HCE)
        THC = (TrackerHitsCollection*)(HCE->GetHC(this->fTrackerCollID));
    else
        return;
    if (THC)
    {
        int n_hit = THC->entries();
        if (n_hit == 0) return;
        this->numberOfEvent ++;
        // whether the interaction point is within the kernel range
        bool validFlag;

        if (this->RecordEventLog)
            G4cout << "/n/n/nEvent ID: " << evt->GetEventID() << G4endl;
        for (int i=0; i<n_hit; i++)
        {
            TrackerHit* aHit = (*THC)[i];
            G4ThreeVector aHitPos = aHit->GetPos();
            // get the position relative to the first interaction point
            aHitPos[2] -= firstInteractionPointZ;
            G4double aHitEdep = aHit->GetEdep();
            size_t idx = this->index(aHitPos, validFlag, i, aHitEdep);
            if (validFlag)
                this->fHitsMap[idx] += aHitEdep;
        }
    }
    else if (this->RecordEventLog)
        G4cout << "No event recorded.\n" << G4endl;
}

void wk::Run::Merge(const G4Run * aRun)
{
    const Run * localRun = static_cast<const Run *>(aRun);
    this->numberOfEvent += localRun->numberOfEvent;
    size_t nElements = this->fHitsMap.size();
    for (int i=0; i<nElements; i++)
        this->fHitsMap[i] += localRun->fHitsMap[i];
}

void wk::Run::logGeom()
{
    // This function logs the geometry parameters
    if (this->pointKernel)
        G4cout << "\n\n\nThis run is in point kernel mode." << G4endl;
    else
        G4cout << "This run is in infinitesimal pencil beam (IPB) mode." << G4endl;
    G4cout << "Below are the geometric parameters" << G4endl;
    G4cout << "The half size of the kernel is (" << this->fSizeX / cm
        << ", " << this->fSizeY / cm << ", " << this->fSizeZ / cm << ") [cm]" << G4endl;
    G4cout << "The dimension of the kernel is  (" << this->fDimX
        << ", " <<  this->fDimY << ", " << this->fDimZ << ")" << G4endl;
    G4cout << "The Z dimension offset is " << this->fPosZ << "\n\n\n" << G4endl;
}