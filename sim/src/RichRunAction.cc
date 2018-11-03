#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ProcessTable.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"

#include "RichRunAction.hh"
#include "RichVAnalysisManager.hh"
#include "RichDetectorConstruction.hh"
#include "RichPrimaryGeneratorAction.hh"
#include "RichSteppingAction.hh"

#include "G4SystemOfUnits.hh"


RichRunAction::RichRunAction(RichVAnalysisManager *an)
{
	analysis = an;
}

RichRunAction::~RichRunAction()
{}

void RichRunAction::BeginOfRunAction(const G4Run * aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

	if (G4VVisManager::GetConcreteInstance()) {
		G4UImanager *UI = G4UImanager::GetUIpointer();
		UI->ApplyCommand("/vis/scene/notifyHandlers");
	}

	RichSteppingAction* steppingAction=(RichSteppingAction*)
		G4RunManager::GetRunManager()->GetUserSteppingAction();

	//Set boundary process pointer to UserSteppingAction
	G4OpBoundaryProcess* opBoundary = (G4OpBoundaryProcess*)
		(G4ProcessTable::GetProcessTable())->FindProcess(G4String("OpBoundary"),
			G4String("opticalphoton"));
	steppingAction->SetBoundaryProcess(opBoundary);

	//Set volume pointers to UserSteppingAction
	RichDetectorConstruction* detector=(RichDetectorConstruction*)
		G4RunManager::GetRunManager()->GetUserDetectorConstruction();
	steppingAction->SetWorldVolume(detector->GetWorldPhysical());
	steppingAction->SetRadiatorVolume(detector->GetRadiatorPhysical());
	steppingAction->SetPMTVolume(detector->GetPMTPhysical());

	RichPrimaryGeneratorAction* primary=(RichPrimaryGeneratorAction*)
		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

#ifdef BGE_ON
	//Set a primary particle as an electron with momentum 1GeV and direction along Z
	//needed for the histrogram initialization
	primary->SetPrimaryName("e-");
	primary->SetTotalMomentum(1*GeV);
	primary->SetTheta(0.);
#endif

    //Set per run primary particle parameters
    primary->InitializeRun();

	analysis->book(aRun);
}

void RichRunAction::EndOfRunAction(const G4Run*)
{
	if (analysis->GetVerboseLevel()>0)
		G4cout << "Run analysis" << G4endl;
	analysis->finish();
}
