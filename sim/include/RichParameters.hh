#ifndef RichParameters_h
# define RichParameters_h 1

# include <vector>

# include "globals.hh"
# include "RichGlobals.hh"

// Limits of Photon Energy  and number of bins for the
// Photon energy range.
static const G4double kPhotonMinEnergy=1.2*CLHEP::eV;
static const G4double kPhotonMaxEnergy=7.4*CLHEP::eV;
static const G4int kNPhotWLbins=100;
extern std::vector<G4double> defPhotMomVector;

static const G4double kPhotonMinWL=kPhotMomWaveConv/(kPhotonMaxEnergy/CLHEP::eV)*CLHEP::nm;
static const G4double kPhotonMaxWL=kPhotMomWaveConv/(kPhotonMinEnergy/CLHEP::RichGeV)*CLHEP::nm;

static const G4double kReferenceWavelength=400*CLHEP::nm;
static const G4double kReferencePhotMom=kPhotMomWaveConv/(kReferenceWavelength/CLHEP::nm)*CLHEP::eV;

// Default radiator geometry & properties
//static const G4double defRadiatorThickness=25*CLHEP::mm;
static const G4double defOverallDimension=100*CLHEP::mm;
//static const G4double defNormThickness1=.055;

static const G4double defAerogelNominalRefIndex=1.05;
static const G4double defAerogelScatLength=5.0*CLHEP::cm;
static const G4double defAerogelSurfTrans=0.95;
static const G4double defProximityDistance=defOverallDimension;
//static const G4int    defNumberOfLayers=4;
//static const G4double defLayer1Thickness=5*CLHEP::mm;

static const G4double kAlpha400=0.438; //coefficient for calculation the aerogel ref.index (cm^3/g):
                                    // n(400)^2 = 1 + alpha_400*rho,  [ NIMA 494(2002)491 ]
static const G4int kMaxNumberOfLayers=50; //largest possible number of layers

static const G4double kDetectorThickness=70*CLHEP::mm; //PMT layer thickness

//PMT default parameters
static const G4double defDetectionEfficiency=1.0;
static const G4double defGeomEfficiency=1.0;

//Paths to input data
#ifdef WORKDIR
static const G4String kWorkDir=WORKDIR;
#else
static const G4String kWorkDir="~/skononov/geant4/farich";
#endif
static const G4String kInputDataDir = kWorkDir+"/input_data/";

static const G4String defAerogelAbsLenDataFile = kInputDataDir+"aer2001_lab.dat";
static const G4String defPMTQEdataFile = kInputDataDir+"hambs.dat";
static const G4String kH2OrefIndDataFile = kInputDataDir+"h2o_ri.dat";
static const G4String kLiFrefIndDataFile = kInputDataDir+"lif_ri.dat";
static const G4String kNaFrefIndDataFile = kInputDataDir+"naf_ri.dat";
static const G4String kCaF2refIndDataFile = kInputDataDir+"caf2_ri.dat";
static const G4String kPMMArefIndDataFile = kInputDataDir+"pmma_ri.dat";

// The origin of the coord system is at the
// upstream end of the World volume.
// The coord system has +z along the beam direction and +y
// going upwards.

// Dimensions of detector parts
static const G4double kWorldHalfX=1200.0*CLHEP::mm;
static const G4double kWorldHalfY=1200.0*CLHEP::mm;
static const G4double kWorldHalfZ=600.0*CLHEP::mm;

static const G4double kRadHalfX=200*CLHEP::mm;
static const G4double kRadHalfY=200*CLHEP::mm;

static const G4double kDetHalfX=1200*CLHEP::mm;
static const G4double kDetHalfY=1200*CLHEP::mm;

// Fused silica refraction index definition
//Malitson, I.H., Journal of the Optical Society of America 55,
//                no. 10 (1965), pp. 1205-1209
static const G4double FusedSilicaDispParameters[6] = {
	0.6961663, 0.0046791, 0.4079426,
	0.0135121, 0.8974794, 97.9340025 };

inline G4double FusedSilicaRefIndex(G4double PhotMom) {
	static const G4double *c = FusedSilicaDispParameters;
	G4double x=pow(kPhotMomWaveConv/(PhotMom/CLHEP::eV)/1e3,2); //wavelength^2 [um^2]
	return sqrt( 1 + x*(c[0]/(x-c[1]) + c[2]/(x-c[3]) + c[4]/(x-c[5])) );
}

static const G4double kReferenceQuartzRefInd=FusedSilicaRefIndex(kReferencePhotMom);

//Parameters of the one-pole Sellmeier formula fit to Novosibirsk aerogel n=1.03 data.
//T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
//Eur. Phys. J. C 52 (2007) 759-764
static const G4double LHCb_a0=0.05639, LHCb_wl0sqr = 83.22*83.22;
static const G4double LHCb_RI2m1ref = LHCb_a0/(1-LHCb_wl0sqr/(400*400));

extern G4double AerogelRefIndexConstant(G4double ri,G4double PhotMom,G4bool getRef=false);

extern G4double AerogelRefIndexQuartz(G4double ri,G4double PhotMom,G4bool getRef=false);

extern G4double AerogelRefIndexLHCb(G4double ri,G4double PhotMom,G4bool getRef=false);

//Aerogel chromaticity models
enum {
	kNoDispersion = 0,
	kQuartzModel  = 1,
	kLHCbModel    = 2
};

extern G4double (*AerogelRefIndex)(G4double, G4double, G4bool getRef=false);

// Function declarations
extern std::vector<G4double> InitializePhotonMomentumVector();

#endif //RichParameters_h

