#include "RichEventAction.hh"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4HCofThisEvent.hh>
#include <G4VHitsCollection.hh>
#include <G4SDManager.hh>
#include <G4VVisManager.hh>
#include "RichVisManager.hh"
#include <G4UImanager.hh>
#include <G4ios.hh>

#include "RichPMT.hh"
#include "RichHit.hh"
#include "RichEventMessenger.hh"
#include "RichTrajectory.hh"

RichEventAction::RichEventAction(RichVAnalysisManager *an) :
	RichCollID(-1),
	forceDrawPhotons(false)
{
	analysis = an;
    messenger = new RichEventMessenger(this);
}

RichEventAction::~RichEventAction() {
	delete messenger;
}

void RichEventAction::BeginOfEventAction(const G4Event *event)
{
	if (RichCollID < 0) {
		RichCollID = G4SDManager::GetSDMpointer()->
			GetCollectionID("RichHitsCollection");
	}

	analysis->BeginOfEventAnalysis(event);
}

void RichEventAction::EndOfEventAction(const G4Event *event)
{
	// first get the trajectories
	G4TrajectoryContainer *trajectoryContainer =
		event->GetTrajectoryContainer();

	G4int n_trajectories = 0;

	G4int verbosity = fpEventManager->GetVerboseLevel();

	if (trajectoryContainer) {
		n_trajectories = trajectoryContainer->entries();
		if (verbosity>0) {
			G4cout << " " << n_trajectories
				   << " tracks are stored in Trajectory Container." << G4endl;
		}
		for (G4int i = 0; i < n_trajectories; i++) {
			RichTrajectory *trj = (RichTrajectory*)(*trajectoryContainer)[i];
			if (trj->GetParticleName()=="opticalphoton") {
				if(forceDrawPhotons) trj->SetDrawTrajectory(true);
			}
			trj->DrawTrajectory(2000);
		}
	}

	// Now get the hits

	if (RichCollID < 0)
		return;

	G4HCofThisEvent *HCE = event->GetHCofThisEvent();
	RichHitsCollection *RHC =
		HCE ? (RichHitsCollection*)HCE->GetHC(RichCollID) : 0;

	G4int n_hit = 0;

	if (RHC) {
		n_hit = RHC->entries();
		for (G4int ih = 0; ih < n_hit; ih++) {
			RichHit *aHit = (*RHC)[ih];
			aHit->SetEvent(event->GetEventID());
			if (aHit->IsHit()) aHit->Draw();
		}
	}

	//Now to Fill the Histograms
	analysis->EndOfEventAnalysis(RHC);
}
