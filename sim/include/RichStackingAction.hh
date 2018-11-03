#ifndef RichStackingAction_h
# define RichStackingAction_h 1

# include <G4UserStackingAction.hh>

class RichStackingAction:public G4UserStackingAction {
public:
	RichStackingAction();
	virtual ~RichStackingAction();

	G4ClassificationOfNewTrack
		ClassifyNewTrack(const G4Track * aTrack);

	void NewStage() {}

	void PrepareNewEvent() {}
};
#endif
