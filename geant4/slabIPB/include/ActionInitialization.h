#ifndef siActionInitialization_h
#define siActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

namespace si
{

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization() = default;
    ~ActionInitialization() override = default;

    // void BuildForMaster() const override;
    void Build() const override;
    void BuildForMaster() const override;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif