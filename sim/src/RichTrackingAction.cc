#include <G4TrackingManager.hh>
#include <G4Track.hh>
#include <G4ParticleTypes.hh>
#include <G4Trajectory.hh>

#include "RichTrackingAction.hh"
#include "RichTrajectory.hh"
#include "RichUserTrackInformation.hh"
#include "RichVAnalysisManager.hh"


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichTrackingAction::PreUserTrackingAction(const G4Track * aTrack)
{
	//Use custom trajectory class
//	fpTrackingManager->SetTrajectory(new RichTrajectory(aTrack));
	fpTrackingManager->SetTrajectory(new G4Trajectory(aTrack));

	//This user track information is only relevant to the photons
	RichUserTrackInformation* trackInfo = new RichUserTrackInformation;
	fpTrackingManager->SetUserTrackInformation(trackInfo);

	if (aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhoton())
	{
  	// particle is optical photon
		const G4VProcess* process = aTrack->GetCreatorProcess();
		if (process)
			if (process->GetProcessName() == "Cerenkov") {
				RichVAnalysisManager::getInstance()->NumCkvPhot++;
				trackInfo->SetLayerNo(aTrack->GetVolume()->GetCopyNo());
                if( aTrack->GetParentID()==1 )
					trackInfo->SetFromPrimary();
			}
	}
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichTrackingAction::PostUserTrackingAction(const G4Track * aTrack)
{
/*	RichTrajectory *trajectory =
		(RichTrajectory *) fpTrackingManager->GimmeTrajectory();

	RichUserTrackInformation *trackInformation =
		(RichUserTrackInformation *) aTrack->GetUserInformation();

	//Lets choose to draw only the photons that hit the PMT
	if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
		if (trackInformation->GetTrackStatus() & trackInformation->hitPMT)
			trajectory->SetDrawTrajectory(true);
	} else //draw all other particles
		trajectory->SetDrawTrajectory(true);
*/
}
