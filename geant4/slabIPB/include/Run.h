#ifndef siRun_h
#define siRun_h 1

#include "G4Run.hh"
#include "G4ThreeVector.hh"

namespace si
{
    class Run : public G4Run
    {
    public:
        Run();
        ~Run() = default;

        virtual void RecordEvent(const G4Event*);
        virtual void Merge(const G4Run*);
        void writeHitsMap(G4String path) const;

        G4float getFullResX() const {return this->fullResX;}
        G4float getFullResY() const {return this->fullResY;}
        G4float getFullResZ() const {return this->fullResZ;}

        G4int getDimX() const {return this->dimX;}
        G4int getDimY() const {return this->dimY;}
        G4int getDimZ() const {return this->dimZ;}

    private:
        std::vector<G4double> fHitsMap;
        // for debug purpose, whether to print the details in RecordEvent
        G4bool recordEventLog;

        // geometric parameters, copied from PhantomDef
        G4float sizeX;
        G4float sizeY;
        G4float sizeZ;
        G4float fullSizeX;
        G4float fullSizeY;
        G4float fullSizeZ;
        G4float resX;
        G4float resY;
        G4float resZ;
        G4float fullResX;
        G4float fullResY;
        G4float fullResZ;
        G4int dimX;
        G4int dimY;
        G4int dimZ;
        std::vector<G4int> HitsCollectionIDs;
        std::vector<G4String> Names;
        inline size_t index(G4ThreeVector position);
    };
}

#endif