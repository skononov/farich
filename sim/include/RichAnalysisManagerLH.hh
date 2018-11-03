#ifndef RichAnalysisManagerLH_h
# define RichAnalysisManagerLH_h 1

# include "globals.hh"
# include "RichVAnalysisManager.hh"
# include "RichParameters.hh"

class G4Run;
class RichHitCollection;

class RichAnalysisManagerLH : public RichVAnalysisManager {
private:
	//private contsructor.
	//use getInstance() to create an object
	RichAnalysisManagerLH();

public:
	virtual ~RichAnalysisManagerLH();

	void book(const G4Run*);
	void finish();

	static RichVAnalysisManager *getInstance() {
		if (instance == 0)
			instance = new RichAnalysisManagerLH;
		return instance;
	}

//	void BeginOfEventAnalysis(const G4Event *);
	void EndOfEventAnalysis(const RichHitsCollection *);

private:
	void LoadBetaData();

//Values per run
	G4double runNumPe;
	G4double runCkvAngle;
	G4double run1peAngleRMS;
	G4double runTrackAngleRMS;
	G4double runRingRadius;
	G4double run1peRadiusRMS;
	G4double runTrackRadiusRMS;
	G4double runBeta;
	G4double run1peBeta;
	G4double run1peBetaRMS;
	G4double runTrackBetaRMS;
	G4int runNumEvWithHit;
	G4int runNumHits;

//Values per event
	G4int layerNumPe[kMaxNumberOfLayers];
	G4double *vRadii, **vLayerRadii; //arrays of event radii
	G4int nRadii; //size of array
};
#endif

