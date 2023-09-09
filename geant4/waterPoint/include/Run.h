#ifndef Run_h
#define Run_h 1

#include <vector>
#include <atomic>

#include "G4Run.hh"
#include "G4THitsMap.hh"
#include "config.h"

namespace wp
{
    class Run : public G4Run
    {
    public:
        Run();
        ~Run();

        virtual void RecordEvent(const G4Event*);
        virtual void Merge(const G4Run*);

        const std::vector<std::vector<std::vector<double>>>&
            GetKernel() const
            {return this->kernel;}
    
    private:
        std::vector<std::tuple<std::string, int, G4THitsMap<G4double>*>> HitsMaps;

        std::vector<std::vector<std::vector<double>>> kernel;
        int PhantomSZ;
        int PhantomBottom;
        // the threshold beyond which the interaction does not count.
        float marginZ;

        int logFrequency;
        int nParticles;
        int PhantomDimXY;
        int PhantomDimZ;
        float resolution;
        float HalfPhantomSizeZ;
    };

    extern std::atomic<int> globalCount;
}

#endif