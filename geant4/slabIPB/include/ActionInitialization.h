#ifndef siActionInitialization_h
#define siActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

namespace si
{
    class ActionInitialization : public G4VUserActionInitialization
    {
    public:
        ActionInitialization();
        ActionInitialization(ActionInitialization& old) = delete;
        ActionInitialization(ActionInitialization&& old) = delete;
        ActionInitialization& operator=(const ActionInitialization& old) = delete;

        ~ActionInitialization();

        virtual void Build() const;
        virtual void BuildForMaster() const;
    };
}

#endif