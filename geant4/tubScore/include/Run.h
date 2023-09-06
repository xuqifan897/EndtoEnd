#ifndef Run_h
#define Run_h 1

#include <vector>

#include "G4Run.hh"
#include "G4THitsMap.hh"

namespace ts
{
    class Run : public G4Run
    {
    public:
        Run();
        ~Run();

        virtual void RecordEvent(const G4Event*);
        virtual void Merge(const G4Run*);
        const std::vector<std::tuple<std::string, int, G4THitsMap<G4double>*>>&
            getHitsMaps() const
            {return this->HitsMaps;}
    
    private:
        // hit collection id
        std::vector<std::tuple<std::string, int, G4THitsMap<G4double>*>> HitsMaps;
    };
}

#endif