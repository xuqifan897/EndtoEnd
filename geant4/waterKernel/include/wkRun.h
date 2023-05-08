#ifndef wkRun_h
#define wkRun_h 1

#include "config.h"
#include "argparse.h"

#include "globals.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include <vector>

#if PHASE > 0
// This class is only needed when we have to accumulate the dose, \\
either for the IPB kernel or the point kernel.

namespace wk
{
    class Run : public G4Run
    {
    public:
        Run();
        virtual ~Run();

        virtual void RecordEvent(const G4Event*);
        virtual void Merge(const G4Run*);

        float getFullSizeX() const {return this->fFullSizeX;}
        float getFullSizeY() const {return this->fFullSizeY;}
        float getFullSizeZ() const {return this->fFullSizeZ;}

        float getFullResX() const {return this->fFullResX;}
        float getFullResY() const {return this->fFullResY;}
        float getFullResZ() const {return this->fFullResZ;}

        int getDimX() const {return this->fDimX;}
        int getDimY() const {return this->fDimY;}
        int getDimZ() const {return this->fDimZ;}

        float getPosZ() const {return this->fPosZ;}
        void writeHitsMap(G4String path) const;
    private:
        std::vector<G4double> fHitsMap;

        // for point kernel, that's the number of valid events, 
        // as some photons pass through the detector without interacting
        // long numberOfEvent;

        // indicates whether it's a point kernel or a IPB kernel
        G4bool pointKernel;
        G4int fTrackerCollID;
        // log the coordinates of individual hits
        G4bool RecordEventLog;

        // half size
        G4float fSizeX;
        G4float fSizeY;
        G4float fSizeZ;

        G4float fFullSizeX;
        G4float fFullSizeY;
        G4float fFullSizeZ;

        // the z position of the particle interaction
        // point w.r.t. this kernel
        G4float fPosZ;
        
        // half size
        G4int fDimX;
        G4int fDimY;
        G4int fDimZ;

        // half size
        G4float fResX;
        G4float fResY;
        G4float fResZ;

        G4float fFullResX;
        G4float fFullResY;
        G4float fFullResZ;

        // the range of the z component of the first interaction point
        // such that the kernel is within the detector. Only used with 
        // point kernel
        G4float fValidRangeZ0;
        G4float fValidRangeZ1;

        // This function is to calculate the index of fHitsMap 
        // from the position
        size_t index(G4ThreeVector position, bool& validFlag, int& hitID, G4double& edep);
        void logGeom();
    };
}

#endif

#endif