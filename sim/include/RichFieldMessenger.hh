#ifndef RichFieldMessenger_h
# define RichFieldMessenger_h 1

# include "globals.hh"
# include "G4UImessenger.hh"

class RichFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class RichMagFieldCommand;

class RichFieldMessenger : public G4UImessenger
{
  public:
	RichFieldMessenger(RichFieldSetup *);
	~RichFieldMessenger();

	void SetNewValue(G4UIcommand *, G4String);
	void SetNewValue(G4UIcommand *, G4int);

  private:

	RichFieldSetup * fEMfieldSetup;

	G4UIdirectory *FieldDir;

	G4UIcmdWithAnInteger *StepperCmd;
	RichMagFieldCommand *MagFieldCmd;
	G4UIcmdWithADoubleAndUnit *MinStepCmd;
	G4UIcmdWithoutParameter *UpdateCmd;
};

class RichMagFieldCommand : public G4UIcommand
{
  public:
	RichMagFieldCommand(const char *theCommandPath,G4UImessenger *theMessenger);

	G4ThreeVector GetFieldVector(const char* paramString);
};

#endif
