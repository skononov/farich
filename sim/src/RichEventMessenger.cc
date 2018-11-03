#include "RichEventMessenger.hh"
#include "RichEventAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichEventMessenger::RichEventMessenger(RichEventAction * evAction) :
	richEvAction(evAction)
{
	forceDrawPhotonsCmd = new G4UIcmdWithABool("/event/forceDrawPhotons", this);
	forceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichEventMessenger::~RichEventMessenger()
{
	delete forceDrawPhotonsCmd;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichEventMessenger::SetNewValue(G4UIcommand * command, G4String newValue)
{
	if (command == forceDrawPhotonsCmd) {
		richEvAction->SetForceDrawPhotons(forceDrawPhotonsCmd->
									  GetNewBoolValue(newValue));
	}
}
