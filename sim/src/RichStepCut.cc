
#include "RichStepCut.hh"

RichStepCut::RichStepCut(const G4String & aName) : G4VDiscreteProcess(aName),
	MaxChargedStep(DBL_MAX)
{
	if(verboseLevel > 0) {
		G4cout << GetProcessName() << " is created " << G4endl;
	}
}

RichStepCut::RichStepCut(RichStepCut & right) : G4VDiscreteProcess(right)
{
}

RichStepCut::~RichStepCut()
{
}

G4double RichStepCut::
PostStepGetPhysicalInteractionLength(const G4Track & aTrack, G4double,
									 G4ForceCondition * condition)
{
	G4double ProposedStep = DBL_MAX;

	if(MaxChargedStep <= 0.)
		return ProposedStep;

	// condition is set to "Not Forced"
	*condition = NotForced;

	if(aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge() == 0)
		return ProposedStep;

	G4VPhysicalVolume* currentVolume=aTrack.GetVolume();

	if( (currentVolume != 0) &&
		(currentVolume->GetName() == "RadiatorPhysical") ||
		(currentVolume->GetMotherLogical() != 0) &&
	    (currentVolume->GetMotherLogical()->GetName() == "RadiatorHolder") )
		ProposedStep = MaxChargedStep;

	return ProposedStep;
}

