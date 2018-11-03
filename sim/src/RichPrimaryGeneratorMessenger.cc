#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

#include "RichPrimaryGeneratorMessenger.hh"
#include "RichPrimaryGeneratorAction.hh"

RichPrimaryGeneratorMessenger::
	RichPrimaryGeneratorMessenger(RichPrimaryGeneratorAction * action) :
	richAction(action)
{
	angleThetaCmd = new G4UIcmdWithADoubleAndUnit("/gun/angleTheta", this);
	angleThetaCmd->SetGuidance("Set polar angle of particle direction.");
	angleThetaCmd->SetParameterName("angle",false);
	angleThetaCmd->SetDefaultUnit("degree");
	angleThetaCmd->SetRange("angle<=180 && angle>=0");
	angleThetaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	anglePhiCmd = new G4UIcmdWithADoubleAndUnit("/gun/anglePhi", this);
	anglePhiCmd->SetGuidance("Set asimutal angle of particle direction.");
	anglePhiCmd->SetParameterName("angle",false);
	anglePhiCmd->SetDefaultUnit("degree");
	anglePhiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	momentumCmd = new G4UIcmdWithADoubleAndUnit("/gun/totalMomentum", this);
	momentumCmd->SetGuidance("Define absolute momentum of the particle.");
	momentumCmd->SetParameterName("momentum",false);
	momentumCmd->SetDefaultUnit("MeV");
	momentumCmd->SetRange("momentum>=0");
	momentumCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	printCmd = new G4UIcommand("/gun/print", this);
	printCmd->SetGuidance("Print particle gun parameters.");
	printCmd->AvailableForStates(G4State_Idle);

	fixVertexCmd = new G4UIcmdWithABool("/gun/fixVertex",this);
	fixVertexCmd->SetGuidance("(Un)Fix primary vertex position.");
	fixVertexCmd->SetParameterName("fix",true);
	fixVertexCmd->SetDefaultValue(true);
	fixVertexCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RichPrimaryGeneratorMessenger::~RichPrimaryGeneratorMessenger()
{
	delete angleThetaCmd;
	delete anglePhiCmd;
	delete momentumCmd;
	delete printCmd;
	delete fixVertexCmd;
}


void RichPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,
												G4String newValue)
{
	if (command == angleThetaCmd)
		richAction->SetTheta(angleThetaCmd->GetNewDoubleValue(newValue));
	else if (command == anglePhiCmd)
		richAction->SetPhi(anglePhiCmd->GetNewDoubleValue(newValue));
	else if (command == momentumCmd)
		richAction->SetTotalMomentum(momentumCmd->GetNewDoubleValue(newValue));
	else if (command == printCmd)
		richAction->Print();
	else if (command == fixVertexCmd)
		richAction->FixVertexPosition(fixVertexCmd->GetNewBoolValue(newValue));
}

