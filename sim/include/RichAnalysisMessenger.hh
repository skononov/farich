#ifndef RichAnalysisMessenger_h
# define RichAnalysisMessenger_h 1

# include "globals.hh"
# include "G4UImessenger.hh"

class RichVAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

class RichAnalysisMessenger : public G4UImessenger {
  public:
	RichAnalysisMessenger(RichVAnalysisManager *);
	~RichAnalysisMessenger();

//	G4String GetCurrentValue(G4UIcommand *);
	void SetNewValue(G4UIcommand *,G4String);

  private:

	//pointer to RichVAnalysisManager
	RichVAnalysisManager *richAnalysis;

	G4UIdirectory 		*RichAnalysisDir;
	G4UIcmdWithAString 	*setInteractiveCommand;
	G4UIcmdWithABool 	*saveHitTreeCommand;
	G4UIcmdWithABool 	*saveEventTreeCommand;
	G4UIcommand 		*runCintCommand;
	G4UIcmdWithAString 	*outputFileCommand;
	G4UIcmdWithAString 	*macroFileCommand;
	G4UIcmdWithAString 	*betaDataFileCommand;
	G4UIcmdWithAnInteger *verboseCommand;
	G4UIcmdWithABool    *weighPixelsCommand;
};
#endif

