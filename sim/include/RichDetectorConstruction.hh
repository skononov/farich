#ifndef RichDetectorConstruction_h
# define RichDetectorConstruction_h 1

# include <G4VUserDetectorConstruction.hh>
# include <G4Material.hh>
# include <G4MaterialTable.hh>
# include <G4MaterialPropertiesTable.hh>
# include <vector>
# include <map>

# include "RichParameters.hh"
# include "RichRunManager.hh"
# include "MLRadiatorDescription.hh"

class RichDetectorMessenger;
class RichPMT;
class RichSpectrum;
class RichAerogelMaterial;
class RichFieldSetup;

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4OpticalSurface;
class G4VisAttributes;

class RichDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	enum {
		fixedNumberOfLayers	=1,
		fixedThickness		=2,
		fixedNT				=4,
		multiRing			=8,
		manual				=16
	}; //MLR parameterization modes
	enum {
		fixedModes			=fixedNumberOfLayers|fixedThickness|fixedNT,
		anyMode				=31
	}; //MLR parameterization masks

	static const char* ModeName(G4int mode);
public:
	RichDetectorConstruction(int verb=0);
	virtual ~RichDetectorConstruction();

	G4VPhysicalVolume* Construct();

	void Update(); //rebuild the detector according to the changes made

	void SetDefaults(); //describe default detector setup with one layer of aerogel

	void Print() const; //print detector description

	void AddLayer(G4String layerSpec,G4double thickness); //add a layer of some material

	void ResetRadiator() { mlrDesc->clear(); Changed(manual); } //reset radiator configuration

// Analytical predictions
	G4double GetNpeRoughly(G4double beta,G4double pathLength);
	G4double GetMeanPhotonMomentum() const { return meanPhotonMomentum; }
private:
	void DefineMaterials();

	void DefineOpPhCathProperties();

	G4LogicalVolume* ConstructRadiator();

	G4VPhysicalVolume* ConstructDetector();

	void SwitchChromaticity();

	void Changed(G4int mask=anyMode) {
		if (mask & mlrMode)
			modified=true;
	}

private:
	RichDetectorMessenger* detectorMessenger;

	RichRunManager* runManager;

	G4VPhysicalVolume* worldPhys;
	G4VPhysicalVolume* radiatorPhys;
	G4VPhysicalVolume* pmtPhys;

	G4LogicalVolume* radiatorLog;

	MLRadiatorDescription* mlrDesc;

	RichFieldSetup *fieldSetup;

	G4MaterialTable aerogelMaterials; //hash table of aerogel materials

	G4int verboseLevel; //verbosity level

	G4bool modified; //is any parameter modified

//Material and geometry parameters
	G4int mlrMode; // MLA parameterization mode

	//parameters independent of MLA calculation mode
	G4double proximityDistance;
	G4double layer1NomRefIndex;
	G4double betaOptimized;

	//radiator distortion parameters
	G4double thicknessDistortion;
	G4double densityDistortion;
	G4double uniformityDistortion;
	G4bool evenDistortion;

	//aerogel optical properties
	G4double aerogelNomScatLength;
	G4double aerogelBoundaryLoss;
	G4String aerogelAbsLenDataFile;
	G4int    chromaticity;

	//fixed thickness mode parameters
	G4double overallDimension;
	G4double normThickness1;

	//fixed number of layers mode parameters
	G4int nLayers;
	G4double layer1Thickness;

	//fixed number of layers & thickness mode
	G4double totalThickness;

	//Multi-ring configuration mode
	G4double series2NomRefIndex;

	//PMT properties
	G4String qeDataFile;
	G4double detectionEfficieny;
	G4double geomEfficiency;
	G4double pixelSize;
	G4double pixelSpacing;
	G4bool   roundPixel;

	//derived quantities
	G4double radiatorThickness;
	G4double radiatorZposition;

	G4double meanPhotonMomentum;

public:
	// Set&get methods
	G4VPhysicalVolume* GetWorldPhysical() { return worldPhys; }
	G4VPhysicalVolume* GetRadiatorPhysical() { return radiatorPhys; }
	G4VPhysicalVolume* GetPMTPhysical() { return pmtPhys; }

	G4int GetVerboseLevel() const { return verboseLevel; }

	G4bool NeedUpdate() const { return modified; }
	G4int GetMLRmode() const { return mlrMode; }

	G4double GetProximityDistance() const { return proximityDistance; }
	G4double GetLayer1NomRefIndex() const { return layer1NomRefIndex; }
	G4double GetBetaOptimized() const { return betaOptimized; }
	G4double GetOverallDimension() const { return overallDimension; }
	G4double GetThickness1Normalized() const { return normThickness1; }
	G4int GetNumberOfLayers() const { return nLayers; }
	G4double GetLayer1Thickness() const { return layer1Thickness; }
	G4double GetSeries2NomRefIndex() const { return series2NomRefIndex; }

	G4double GetAerogelScatLength() const { return aerogelNomScatLength; }
	G4double GetAerogelBoundaryLoss() const { return aerogelBoundaryLoss; }
	G4String GetAerogelAbsorptionFile() const { return aerogelAbsLenDataFile; }
	G4bool IsChromaticityOn() const { return chromaticity!=0; }
	G4int Chromaticity() const { return chromaticity; }

	G4double GetDetectionEfficiency() const { return detectionEfficieny; }
	G4String GetQEfile() const { return qeDataFile; }
	G4double GetGeomEfficiency() const { return geomEfficiency; }
	G4double GetPixelSize() const { return pixelSize; }
	G4double GetPixelSpacing() const { return pixelSpacing; }

	G4double GetRadiatorThickness() const { return radiatorThickness; }
	G4double GetRadiatorZposition() const { return radiatorZposition; }
	G4double GetRadiatorOutSurfZ() const { return radiatorZposition+radiatorThickness/2; }
	G4double GetRadiatorInSurfZ() const { return radiatorZposition-radiatorThickness/2; }
	G4double GetDetectorWindowZ() const { return GetRadiatorOutSurfZ()+proximityDistance; }

	RichFieldSetup *GetFieldSetup() const { return fieldSetup; }

	MLRadiatorDescription* GetMLRadiatorDescription() const { return mlrDesc; }

	void SetVerboseLevel(G4int level) { verboseLevel=level; }

	void SetMLRmode(G4int mode) {
		if (mode!=mlrMode) Changed();
		mlrMode=mode;
	}
	void SetProximityDistance(G4double d) {
		proximityDistance=d; Changed(~fixedNT&anyMode);
	}
	void SetLayer1NomRefIndex(G4double ri) {
		layer1NomRefIndex=ri; Changed();
	}
	void SetBetaOptimized(G4double b) {
		betaOptimized=b;
		Changed(~multiRing&anyMode);
	}
	void SetOverallDimension(G4double d) {
		overallDimension=d; Changed(fixedThickness|fixedNT|multiRing);
	}
	void SetThickness1Normalized(G4double t) {
		normThickness1=t; Changed(fixedThickness);
	}

	void SetNumberOfLayers(G4int n) {
		nLayers = n>kMaxNumberOfLayers?kMaxNumberOfLayers:n;
		Changed(fixedNumberOfLayers|fixedNT|multiRing);
	}
	void SetLayer1Thickness(G4double t) {
		layer1Thickness=t;
		Changed(fixedNumberOfLayers);
	}
	void SetTotalThickness(G4double t) {
		totalThickness=t;
		Changed(fixedNT);
	}
	void SetSeries2NomRefIndex(G4double ri) {
		series2NomRefIndex = ri;
		Changed(multiRing);
	}

	void SetThicknessDistortion(G4double eps) {
		thicknessDistortion=eps;
		Changed();
	}
	void SetDensityDistortion(G4double eps) {
		densityDistortion=eps;
		Changed();
	}
	void SetUniformityDistortion(G4double eps) {
		uniformityDistortion=eps;
		Changed();
	}
	void SetEvenDistortion(G4bool on=true) {
		evenDistortion=on;
		Changed();
	}

	void SetAerogelScatLength(G4double  lsc) {
		aerogelNomScatLength=lsc; Changed();
	}
	void SetAerogelBoundaryLoss(G4double bl) {
		aerogelBoundaryLoss=bl; Changed();
	}
	void SetAerogelAbsorptionFile(const G4String& fn) {
		aerogelAbsLenDataFile=fn; Changed();
	}
	void SetDetectionEfficiency(G4double eff) {
		detectionEfficieny=eff; Changed();
	}
	void SetGeomEfficiency(G4double eff) { geomEfficiency=eff; Changed(); }
	void SetPixelSize(G4double size) { pixelSize=size; Changed(); }
	void SetQEfile(const G4String& fn) { qeDataFile=fn; Changed(); }
	void SetChromaticity(G4int model) { chromaticity=model; Changed(); }

private:
//Materials and surfaces
	G4MaterialTable RichUsedMaterials;
	G4OpticalSurface* RichOpPhotocathodeSurface;

//Properties tables
	G4MaterialPropertiesTable* RichAirMPT;
	G4MaterialPropertiesTable* RichQuartzMPT;
	G4MaterialPropertiesTable* RichH2OMPT;
	G4MaterialPropertiesTable* RichLiFMPT;
	G4MaterialPropertiesTable* RichNaFMPT;
	G4MaterialPropertiesTable* RichCaF2MPT;
	G4MaterialPropertiesTable* RichPMMAMPT;
	RichAerogelMaterial*       RichAerogelMat;
	G4MaterialPropertiesTable* RichOpPhCathSurfacePT;

	RichPMT* pmtSD;

    typedef std::map<G4Material*,G4MaterialPropertiesTable*> MMPT_t;
	MMPT_t DefaultMPTs;

//Visual attributes
	G4VisAttributes* WorldVisAttr;
	G4VisAttributes* AerogelVisAttr;
	G4VisAttributes* PMTVisAttr;
};

inline const char* RichDetectorConstruction::ModeName(G4int mode)
{
	if (mode==fixedNumberOfLayers) return "fixed_layers_number";
	if (mode==fixedThickness)      return "fixed_thickness";
	if (mode==fixedNT)             return "fixed_all";
	if (mode==multiRing)           return "multiring";
	if (mode==manual)              return "manual";
	return "unknown";
}

#endif
