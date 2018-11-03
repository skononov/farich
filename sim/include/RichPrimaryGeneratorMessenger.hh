#ifndef RichPrimaryGeneratorMessenger_h
# define RichPrimaryGeneratorMessenger_h 1

# include "G4UImessenger.hh"
# include "globals.hh"

class RichPrimaryGeneratorAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

class RichPrimaryGeneratorMessenger : public G4UImessenger {
public:
	RichPrimaryGeneratorMessenger(RichPrimaryGeneratorAction *);
	~RichPrimaryGeneratorMessenger();

	void SetNewValue(G4UIcommand *, G4String);

private:
	RichPrimaryGeneratorAction * richAction;
	G4UIcmdWithADoubleAndUnit *angleThetaCmd;
	G4UIcmdWithADoubleAndUnit *anglePhiCmd;
	G4UIcmdWithADoubleAndUnit *momentumCmd;
	G4UIcmdWithABool          *fixVertexCmd;
	G4UIcommand               *printCmd;
};

#endif
