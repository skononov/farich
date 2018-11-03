#ifndef RichAnalysisManagerS_h
# define RichAnalysisManagerS_h 1

# include "globals.hh"

# include "RichParameters.hh"
# include "RichVAnalysisManager.hh"

class TF1;
class G4Run;

class RichAnalysisManagerS : public RichVAnalysisManager {
private:
	//private contsructor.
	//use getInstance() to create an object
	RichAnalysisManagerS();

public:
	virtual ~RichAnalysisManagerS();

	void book(const G4Run*);
	void finish();

	static RichVAnalysisManager *getInstance() {
		if (instance == 0)
			instance = new RichAnalysisManagerS;
		return instance;
	}

//	void BeginOfEventAnalysis(const G4Event *);
	void EndOfEventAnalysis(const RichHitsCollection *);

private:
	void LoadBetaData();

//radius to beta transformation
	TF1 *fBeta;
	TF1 *fSigmaBeta;
	TF1 *fRing;

//Values per run
	G4double runNumPe;
	G4double runCkvAngle;
	G4double run1peAngle;
	G4double run1peAngleRMS;
	G4double runTrackAngleRMS;
	G4double runRingRadius;
	G4double run1peRadius;
	G4double run1peRadiusRMS;
	G4double runTrackRadiusRMS;
	G4double runBeta;
	G4double run1peBetaRMS;
	G4double runTrackBetaRMS;
	G4int runNumEvWithHit;
	G4int runNumHits;

//Values per event
	G4int layerNumPe[kMaxNumberOfLayers];
	G4double layerBeta[kMaxNumberOfLayers];
};
#endif

