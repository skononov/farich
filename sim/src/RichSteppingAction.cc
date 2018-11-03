#include <G4SDManager.hh>
#include <G4EventManager.hh>
#include <G4ProcessManager.hh>
#include <G4Track.hh>
#include <G4Step.hh>
#include <G4Event.hh>
#include <G4StepPoint.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4OpBoundaryProcess.hh>

#include "RichSteppingAction.hh"
#include "RichPMT.hh"
#include "RichUserTrackInformation.hh"
#include "RichVAnalysisManager.hh"

#include "G4SystemOfUnits.hh"


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichSteppingAction::UserSteppingAction(const G4Step * theStep)
{
	G4Track *theTrack = theStep->GetTrack();

	G4ParticleDefinition *particleType = theTrack->GetDefinition();

	RichUserTrackInformation *trackInfo
		= (RichUserTrackInformation *) theTrack->GetUserInformation();

	RichVAnalysisManager* analysis = RichVAnalysisManager::getInstance();

	G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();

	G4StepPoint *thePostPoint = theStep->GetPostStepPoint();
	G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();

	if (!thePostPV)	return;		//out of world

	if (particleType!=G4OpticalPhoton::OpticalPhoton() && theTrack->GetKineticEnergy()<0.001*MeV) {
		theTrack->SetTrackStatus(fStopAndKill);
		return;
	}

	G4bool isAtBoundary=(thePostPoint->GetStepStatus()==fGeomBoundary);

#ifndef BGE_ON
	if (theTrack->GetTrackID()==1 && isAtBoundary)
	{
		if (thePrePV == worldVolume) {
			G4VPhysicalVolume *thePostPVmother =
				thePostPoint->GetTouchable()->GetVolume(1);
			if (thePostPV == radiatorVolume || thePostPVmother == radiatorVolume)
				analysis->inPrimaryPosition=theTrack->GetPosition();
		} else if (thePostPV == worldVolume) {
			G4VPhysicalVolume *thePrePVmother =
				thePrePoint->GetTouchable()->GetVolume(1);
			if (thePrePV == radiatorVolume || thePrePVmother == radiatorVolume)
				analysis->outPrimaryPosition=theTrack->GetPosition();
		}
		return;
	}
#else
	if (theTrack->GetTrackID()==1) return;
#endif

	G4OpBoundaryProcessStatus boundaryStatus = Undefined;
	static RichPMT *pmt = NULL;

	pmt = (RichPMT*)G4SDManager::GetSDMpointer()->
		FindSensitiveDetector(G4String("/RichPMT"));

	if (!opBoundary) return;

	if (particleType == G4OpticalPhoton::OpticalPhoton()) {
		//Optical photon only

		//Was the photon absorbed by the absorption process
		const G4VProcess* process = thePostPoint->GetProcessDefinedStep();
		G4String pname = process->GetProcessName();

		if (pname == "OpAbsorption")
			trackInfo->AddTrackStatusFlag(trackInfo->absorbed);
		else if (pname == "OpRayleigh") {
			trackInfo->IncScatterCount();
#ifndef BACKGROUND_ON
            //do not track scattered photon further
			theTrack->SetTrackStatus(fStopAndKill);
#endif
		}

		//Check to see if the particle was actually at a boundary
		//Otherwise the boundary status may not be valid
		//Prior to Geant4.6.0-p1 this would not have been enough to check
		if (isAtBoundary) {
//			G4String PreVolName = thePrePV->GetName();
//			G4String PostVolName = thePostPV->GetName();

			G4double incangle=0;
			if ( thePostPV == pmtVolume ) {
				analysis->NumPhotAtPMT++;
				incangle=(thePrePoint->GetMomentumDirection()).theta()/deg;
				(analysis->hNumPcIncidenceVsAngle)->Fill(incangle);
			}

			boundaryStatus = opBoundary->GetStatus();
			switch (boundaryStatus) {
				case Absorption:
					trackInfo->AddTrackStatusFlag(trackInfo->boundaryAbsorbed);
					break;
				case Detection:
					//Note, this assumes that the volume causing detection
					//is the photocathode because it is the only one with
					//non-zero efficiency
				{
					//Trigger sensitive detector manually since photon is
					//absorbed but status was Detection
//					trackInfo->AddTrackStatusFlag(trackInfo->hitPMT);
                    if (pmt) pmt->ProcessHits_constStep(theStep);
					break;
				}
				case TotalInternalReflection:
				case SpikeReflection:
//					trackInfo->IncReflectionCount();
					theTrack->SetTrackStatus(fStopAndKill);
                    break;
				case FresnelReflection:
					if (thePostPV == pmtVolume)
						(analysis->hPcReflectionVsAngle)->Fill(incangle);
#ifndef BACKGROUND_ON
					theTrack->SetTrackStatus(fStopAndKill);
#endif
                    break;
				default:
					break;
			} //switch(boundaryStatus)
		} //if on boundary
	} //if optical photon
}
