#ifndef RichRunAction_h
# define RichRunAction_h 1

# include "globals.hh"
# include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class RichVAnalysisManager;

class RichRunAction : public G4UserRunAction {
public:
	RichRunAction(RichVAnalysisManager *);
	~RichRunAction();

public:
	void BeginOfRunAction(const G4Run * aRun);
	void EndOfRunAction(const G4Run * aRun);

private:
	RichVAnalysisManager *analysis;
};

#endif
