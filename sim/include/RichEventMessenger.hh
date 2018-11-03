#ifndef RichEventMessenger_hh
# define RichEventMessenger_hh 1

# include "G4UImessenger.hh"
# include "globals.hh"

class RichEventAction;
class G4UIcmdWithABool;

class RichEventMessenger : public G4UImessenger {
public:
	RichEventMessenger(RichEventAction *);
	~RichEventMessenger();

	void SetNewValue(G4UIcommand *, G4String);

private:
	RichEventAction *richEvAction;
	G4UIcmdWithABool *forceDrawPhotonsCmd;
};

#endif
