#ifndef RichSteppingAction_hh
# define RichSteppingAction_hh 1

# include "globals.hh"
# include <G4UserSteppingAction.hh>

class G4OpBoundaryProcess;
class G4VPhysicalVolume;

class RichSteppingAction : public G4UserSteppingAction {
public:
	RichSteppingAction() :
		opBoundary(0), worldVolume(0), radiatorVolume(0), pmtVolume(0)
	{}
	~RichSteppingAction() {}

	void UserSteppingAction(const G4Step *);

private:
	G4OpBoundaryProcess* opBoundary;
	G4VPhysicalVolume* worldVolume;
	G4VPhysicalVolume* radiatorVolume;
	G4VPhysicalVolume* pmtVolume;

public:
	void SetBoundaryProcess(G4OpBoundaryProcess* process)
	{ opBoundary=process; }
	void SetWorldVolume(G4VPhysicalVolume* volume)
	{ worldVolume=volume; }
	void SetRadiatorVolume(G4VPhysicalVolume* volume)
	{ radiatorVolume=volume; }
	void SetPMTVolume(G4VPhysicalVolume* volume)
	{ pmtVolume=volume; }
};

#endif
