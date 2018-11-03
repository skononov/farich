#include "RichStackingAction.hh"

#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>

RichStackingAction::RichStackingAction()
{}

RichStackingAction::~RichStackingAction()
{}

G4ClassificationOfNewTrack
RichStackingAction::ClassifyNewTrack(const G4Track *aTrack)
{
	if (aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
	{
  	// particle is optical photon
		return fUrgent;
	}
    return fWaiting;
}

