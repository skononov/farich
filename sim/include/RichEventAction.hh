#ifndef RichEventAction_h
# define RichEventAction_h 1

# include <G4UserEventAction.hh>
# include <G4ThreeVector.hh>
# include "RichVAnalysisManager.hh"

class G4Event;
class RichEventMessenger;

class RichEventAction : public G4UserEventAction {
public:
	RichEventAction(RichVAnalysisManager *);
	virtual ~RichEventAction();

public:
	void BeginOfEventAction(const G4Event *);
	void EndOfEventAction(const G4Event *);
	G4int GetRichCollID() {	return RichCollID; }
	void SetForceDrawPhotons(G4bool b=true) { forceDrawPhotons=b; }

private:
	G4int RichCollID;
	RichVAnalysisManager *analysis;
	G4bool forceDrawPhotons;
	RichEventMessenger *messenger;
};
#endif
