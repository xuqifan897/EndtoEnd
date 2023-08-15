#include "StepMax.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ts::StepMax::StepMax(const G4String& processName)
: G4VDiscreteProcess(processName), fMaxChargedStep(DBL_MAX)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ts::StepMax::~StepMax() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ts::StepMax::IsApplicable(const G4ParticleDefinition& particle)
{
    return (particle.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ts::StepMax::SetMaxStep(G4double step) {fMaxChargedStep = step;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ts::StepMax::PostStepGetPhysicalInteractionLength(const G4Track&,
                                                       G4double,
                                                       G4ForceCondition* condition )
{
    // condition is set to "Not Forced"
    *condition = NotForced;
    
    G4double ProposedStep = DBL_MAX;
    
    if(fMaxChargedStep > 0.) ProposedStep = fMaxChargedStep;
    
    return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* ts::StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
    // do nothing
    aParticleChange.Initialize(aTrack);
    return &aParticleChange;
}