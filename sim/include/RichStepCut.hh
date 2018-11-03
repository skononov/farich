#ifndef RichStepCut_h
# define RichStepCut_h 1

# include "G4ios.hh"
# include "globals.hh"
# include "G4VDiscreteProcess.hh"
# include "G4Track.hh"
# include "G4Step.hh"
# include "G4VParticleChange.hh"
# include "G4EnergyLossTables.hh"
# include "G4UserLimits.hh"

class RichStepCut : public G4VDiscreteProcess
{
public:

	RichStepCut(const G4String & processName = "UserStepCut");
	RichStepCut(RichStepCut &);

	~RichStepCut();

	G4bool IsApplicable(const G4ParticleDefinition& p);

	G4double PostStepGetPhysicalInteractionLength(const G4Track & track,
												  G4double previousStepSize,
												  G4ForceCondition * condition);

	G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &);

	void SetMaxStep(G4double);

protected:
	// it is not needed here !
	G4double GetMeanFreePath(const G4Track & aTrack,
							   G4double previousStepSize,
							   G4ForceCondition * condition);

private:
	// hide assignment operator as private
	RichStepCut & operator=(const RichStepCut & right);

private:
	G4double MaxChargedStep;
};

// inlined function members implementation


inline G4VParticleChange *RichStepCut::PostStepDoIt(const G4Track & aTrack,
												   const G4Step &)
{
	// do nothing
	aParticleChange.Initialize(aTrack);
	return &aParticleChange;
}

inline G4double RichStepCut::GetMeanFreePath(const G4Track &,
											G4double, G4ForceCondition *)
{
	return 0.;
}

inline G4bool RichStepCut::IsApplicable(const G4ParticleDefinition& p)
{
	if( p.GetPDGCharge() != 0.0 ) return true;
	return false;
}

inline void RichStepCut::SetMaxStep(G4double step)
{
	MaxChargedStep = step;
}

#endif
