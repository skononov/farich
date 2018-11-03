#ifndef RichTrackingAction_h
# define RichTrackingAction_h 1

# include <G4UserTrackingAction.hh>

class RichTrackingAction : public G4UserTrackingAction {
public:
	RichTrackingAction() {}
	~RichTrackingAction() {}

	void PreUserTrackingAction(const G4Track *);
	void PostUserTrackingAction(const G4Track *);

};

#endif
