#ifndef RichVAnalysisManager_h
# define RichVAnalysisManager_h 1

# include "globals.hh"
# include <cmath>
# include <G4ThreeVector.hh>
# include "Rtypes.h"
# include "TRint.h"
# include "TH2.h"
# include "RichHit.hh"
# include "RichParameters.hh"
# include "RichHitInfo.hh"

class G4Step;
class G4Run;
class G4Event;
class G4Step;
class G4Timer;
class G4SteppingManager;
class RichAnalysisMessenger;
class RichDetectorConstruction;
class RichPrimaryGeneratorAction;
class RichRootEvent;

class TROOT;
class TFile;
class TTree;
class TProfile;

class RichVAnalysisManager {
protected:
	//protected contsructor.
	RichVAnalysisManager();

public:
	virtual ~RichVAnalysisManager();

	virtual void book(const G4Run*);
	virtual void finish();

	static RichVAnalysisManager *getInstance() {
		return instance;
	}

	virtual void BeginOfEventAnalysis(const G4Event *);
	virtual void EndOfEventAnalysis(const RichHitsCollection *)=0;

	void RunCint() { application->Run(kTRUE); }

protected:
	G4bool PropagateParticle(G4double z0, G4ThreeVector& endPosition, G4ThreeVector& endDirection);

	inline G4double tanrefangle(G4double tan_i,G4double n_i,G4double n_k);

	void find_thetac(const G4ThreeVector& hitpos,G4int iLayer,G4double& thetac,G4double& alpha);

	virtual void LoadBetaData()=0;

	static RichDetectorConstruction* detector;
	static RichPrimaryGeneratorAction* primary;

	static RichVAnalysisManager *instance;

	std::ostream ancout;

	RichAnalysisMessenger *analysisMessenger;

	TRint *application;

	G4Run *currentRun;

	TFile *outputFile; //root file to store histograms

	G4Timer *iTimer;

//Analysis configuration parameters
	G4String outputFileName;    //root file name, rich.root by default
	G4String macroFileName;     //macro file name executed at the end of each run, run.C by default
	G4String betaDataFileName;  //beta calibration file name, no by default
	G4bool interactive;         //if analysis interactive, true by default
	G4int verboseLevel;         //verbosity level, 1 by default
	G4bool saveHitTree;         //hit tree storing flag, true by default
	G4bool saveEventTree;       //event tree storing flag, false by default
	G4bool weighPixels; 		//weigh pixels by number of photoelectrons, false by default

public:
	void SetOutputFileName(const G4String& newName) { outputFileName = newName; }
	void SetMacroFileName(const G4String& newName) { macroFileName = newName; }
	void SetBetaDataFileName(const G4String& newName)
	{ betaDataFileName = newName; }
	void SetInteractive(G4bool b=true) { interactive=b; }
	void SetVerboseLevel(G4int vl) { verboseLevel=vl; }
	void SetHitTreeStoring(G4bool b=true) { saveHitTree=b; }
	void SetEventTreeStoring(G4bool b=true) { saveEventTree=b; }
	void SetWeighPixels(G4bool b=true) { weighPixels=b; }

	const G4String& GetOutputFileName() const { return outputFileName; }
	const G4String& GetMacroFileName() const { return macroFileName; }
	const G4String& SetBetaDataFileName() const { return betaDataFileName; }
	G4bool IsInteractive() const { return interactive; }
	G4int GetVerboseLevel() const { return verboseLevel; }
	G4bool GetHitTreeStoring() const { return saveHitTree; }
	G4bool GetEventTreeStoring() const { return saveEventTree; }
	G4bool GetWeighPixels() const { return weighPixels; }

public:
//Values per event
	G4int NumCkvPhot;
	G4int NumPhotoelectrons;
	G4int NumPhotAtPMT;

//Points where primary crosses the radiator faces. Set by RichSteppingAction for each event.
	G4ThreeVector inPrimaryPosition, outPrimaryPosition;

//Magnetic field value
	G4ThreeVector magField;

//Distributions for a run
	TH1F *hNumPcIncidenceVsAngle;
	TH1F *hPcReflectionVsAngle;

protected:
//Used by derived objects
	G4double expectedNpe;

//Histograms per run
	TH1F *hNumCkvPhot;
	TH1F *hNumPhotAtPMT;
	TH1F *hNumPhotoelectrons;
	TH1F *hPeWL;

	TH2F *hXYhits;
	TH1F *h1peBeta; //beta distribution per track
	TH1F *hCkvAngle; //Cherenkov angle distrubution per photon
	TH1F *hAlphaAngle; //Azimuthal angle distrubution per photon
    TProfile *pCkvAngVsAlpha; //Cherenkov angle as functions of alpha angle
	TH1F *h1peRadius; //Cherenkov radius distrubution per photon
	TH1F *hLayerNumPe[kMaxNumberOfLayers]; //Number of photoelectrons per layer
	TH1F *hLayerCkvAngle[kMaxNumberOfLayers]; //Cherenkov angle distrubution per photon per layer
	TH1F *hLayer1peRadius[kMaxNumberOfLayers]; //Cherenkov radius distrubution per photon per layer
	TH1F *hLayer1peBeta[kMaxNumberOfLayers]; //Measured beta distrubution per photon per layer
	TH1F *hLayerTrackBeta[kMaxNumberOfLayers]; //Measured beta distrubution per track per layer
	TH1F *hTrackCkvAngle; //Cherenkov angle distrubution per track
	TH1F *hTrackRadius; //Cherenkov radius distrubution per track
	TH1F *hTrackBeta; //Measured beta distrubution per track
	TH1F *hNpePerPixel; //Distribution of number of photoelectrons per pixel

	TTree *tHits; //Stores information per hit
	TTree *tEvent; //Stores events

	struct HitInfo hitinfo; //information per hit

	RichRootEvent* rootEvent; //event to store in the tree

//Auxilary parameters
	G4bool betaDataExist; //if beta data successfully read
	G4double beta;
	G4double momentum; //primary particle momentum
	G4double charge; //primary particle charge
	G4double layerIndices[kMaxNumberOfLayers]; //Radiator layer indices
	G4double layerThicknesses[kMaxNumberOfLayers]; //Radiator layer thicknesses
	G4double layerAngles[kMaxNumberOfLayers]; //Cherenkov angles in radiator layers
	G4double proximityDistance; //proximity distance
	G4double incAngle; //particle incident angle counted from normal, rad
	G4int nLayers; //number of aerogel layers
	G4int nRings; //number of rings
	G4int mlrMode; //MLR parameterization mode
	G4ThreeVector primVertexPosition; //primary particle vertex position
	G4ThreeVector primDirection; //primary particle direction
	G4ThreeVector meanEmissionPoint; //Cherenkov photon mean vertex position
	G4ThreeVector meanPrimDirection; //Mean direction of the primary particle in the radiator
	G4ThreeVector layerVertexPositions[kMaxNumberOfLayers]; //Positions of mean photon vertex in radiator layers
	G4ThreeVector eventVertexShift; //shift of the primary vertex in an event
	G4ThreeVector imageCenter; //particle hit at detection plane
};

inline G4double RichVAnalysisManager::tanrefangle(G4double tan_i,G4double n_i,G4double n_k)
{
	G4double n=n_i/n_k;
	G4double denom2=1+(1-n*n)*tan_i*tan_i;
	if (denom2<=0) return NAN;

	return n*tan_i/sqrt(denom2);
}

#endif

