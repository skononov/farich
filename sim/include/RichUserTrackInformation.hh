#ifndef RichUserTrackInformation_hh
# define RichUserTrackInformation_hh 1

#include <G4VUserTrackInformation.hh>
#include "globals.hh"

/*RichTrackStatus:
  active: still being tracked
  hitPMT: stopped by being detected in a PMT
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  scattered: undergone Rayleigh scattering process
  reflected: undergone a Fresnel reflection
  inactive: track is stopped for some reason
  -This is the sum of all stopped flags so can be used to remove stopped flags
*/

class RichUserTrackInformation:public G4VUserTrackInformation {
public:
	enum RichTrackStatus {
		active = 1,
		hitPMT = 2,
		absorbed = 4,
		boundaryAbsorbed = 8,
		inactive = 14
	};

public:
	RichUserTrackInformation() :
		status(active), nscatter(0), nrefl(0), fromPrimary(false)
	{}
	~RichUserTrackInformation()
	{}

	//Sets the track status to s (does not check validity of flags)
	void SetTrackStatusFlags(int s) { status = s; }
	//Does a smart add of track status flags (disabling old flags that conflict)
	//If there's conflicts with itself it will not be detected
	void AddTrackStatusFlag(int s)
	{
		if (s & active)				//track is now active
			status &= ~inactive;	//remove any flags indicating it is inactive
		else if (s & inactive)		//track is now inactive
			status &= ~active;		//remove any flags indicating it is active
		status |= s;				//add new flags
	}

	int GetTrackStatus() const { return status;	}

	void IncScatterCount() { nscatter++; }
	G4int GetScatterCount() const { return nscatter; }

	void IncReflectionCount() { nrefl++; }
	G4int GetReflectionCount() const { return nrefl; }

	void SetLayerNo(G4int n) { layerNo = n; }
	G4int GetLayerNo() const { return layerNo; }

	void SetFromPrimary(G4bool p=true) { fromPrimary=p; }
	G4bool IsFromPrimary() const { return fromPrimary; }

	void Print() const {}

private:
	int status;
	G4int nscatter, nrefl;
	G4int layerNo;
	G4bool fromPrimary;
};

#endif
