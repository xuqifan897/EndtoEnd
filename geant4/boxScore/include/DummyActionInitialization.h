#ifndef DummyActionInitialization_h
#define DummyActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

namespace bs
{
    class DummyActionInitialization : public G4VUserActionInitialization
    {
    public:
        DummyActionInitialization() = default;
        ~DummyActionInitialization() override = default;

        void BuildForMaster() const override;
        void Build() const override;
    };
}

#endif