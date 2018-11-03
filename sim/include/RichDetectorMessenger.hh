#ifndef RichDetectorMessenger_h
#define RichDetectorMessenger_h 1

#include <G4UImessenger.hh>
#include <G4UIcommand.hh>
#include "globals.hh"

class RichDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

class AddLayerCommand; //see below

class RichDetectorMessenger : public G4UImessenger
{
public:
	RichDetectorMessenger(RichDetectorConstruction*);
	~RichDetectorMessenger();

	void SetNewValue(G4UIcommand*, G4String);

private:
	RichDetectorConstruction* 	richDetector;

// Global setup
	G4UIdirectory*              richDir;              // detector setup directory
	G4UIcmdWithAnInteger*       verboseCmd;           // set verbosity
	G4UIcmdWithoutParameter*	defaultsCmd;          // reset to defaults
	G4UIcmdWithoutParameter*	updateCmd;            // build modified detector
	G4UIcmdWithoutParameter*	printCmd;             // print detector's setup
	G4UIcmdWithAString*			selectModeCmd;        // radiator constructing mode selector
	/* There are four predefined radiator constructing modes:
       "r" - Radiator is to be defined manually with 'addAerogel' and/or 'addMaterial' commands.
	   "n" - Single ring with fixed number of layers. The radiator configuration is defined by
	         proximity distance along with number of layers.
	   "t" - Single ring with fixed total thickness of radiator. Thickness is determined from
	         difference of overall detector dimension and proximity distance.
	   "f" - Single ring with fixed both total radiator thickness and number of layers. Thickness
			 of radiator is given directly by the 'totalThickness' command.
	   "m" - Multi-ring mode with two series of layers with different indices and the same
	         thickness. One need to specify proximity distance, overall dimension and number
			 of layers.

	   Refractive index of the first layer is to be set in all modes with layer1RefIndCmd.
	   In n-mode one need to give the absolute thickness of the first layer by layer1ThicknessCmd,
	   while in t-mode - ratio of first layer thickness to proximity distance by normThickness1Cmd.
	   Constructing radiator in single ring modes (n and t) can be performed for given beta
	   by betaOptimizedCmd.
	   In m-mode second series index is set by series2RefIndCmd.
	*/

	G4UIcmdWithADoubleAndUnit*	proximityDistanceCmd; // set span distance between radiator face and PMT
	G4UIcmdWithADoubleAndUnit*	overallDimensionCmd;  // set overall size used by radiator and span

// Radiator setup
	G4UIdirectory*              radiatorDir;          // radiator setup directory

	// define radiator layers manually
	AddLayerCommand*			addLayerCmd;          // add a layer of radiator
	G4UIcmdWithoutParameter*	resetCmd;             // reset previously defined layers

	// complex radiator constructing parameters valid in n,t,f,m modes
	G4UIcmdWithAnInteger*		nLayersCmd;           // set number of layers  (n,f,m)
	G4UIcmdWithADouble*			layer1RefIndCmd;      // set index of the first layer (all)
	G4UIcmdWithADoubleAndUnit*	layer1ThicknessCmd;   // set thickness of the first layer (n)
	G4UIcmdWithADouble*			normThickness1Cmd;    // thickness of the first layer
                                                      //  in fractions of proximity distance (t)
	G4UIcmdWithADoubleAndUnit*	totalThicknessCmd;    // set total thickness of radiator (f)

	G4UIcmdWithADouble*			betaOptimizedCmd;     // optimise radiator focusing for given beta (n,t,f)
	G4UIcmdWithADouble*			series2RefIndCmd;     // set index of second series layers (m)

	// distortion parameters
	G4UIcmdWithADouble*			distortThicknessCmd;  // distort thickness of every layer (except r)
	G4UIcmdWithADouble*			distortIndexCmd;      // distort ref. index of every layer (except r)
	G4UIcmdWithADouble*			distortUniformityCmd; // distort ref. index uniformity of every layer (except r)
	G4UIcmdWithABool*			evenDistortionCmd;    // set distortion to be even or odd (except r)

    // common optic properties of aerogel and other materials
	G4UIcmdWithADoubleAndUnit*	scatteringLengthCmd;  // set Rayleigh scattering length of aerogel at 400nm
	G4UIcmdWithAString*         absLengthDataFileCmd; // set name of file with absorption length data
	G4UIcmdWithADouble*			boundaryLossCmd;      // set boundary loss fraction
	G4UIcmdWithAString*			chromaticityCmd;      // choice of chromaticity model

// Photodetector setup
	G4UIdirectory*              pmtDir;               // PMT setup directory
	G4UIcmdWithAString*         qeDataFileCmd;        // set file with QE data of PMT
	G4UIcmdWithADouble*			detectionEffCmd;      // set detection efficiency of the detector
	G4UIcmdWithADouble*			geomEffCmd;           // set geometrical efficiency of the detector (pixel density)
	G4UIcmdWithADoubleAndUnit*	pixelSizeCmd;         // set pixel size
};

// My command for layer insertion. Takes three arguments:
//  index/material_name, thickness, unit of thickness.
// Type of the first argument is not predefined in this class. Two type of objects can define two
// different commands: first argument may be a number or a string. In the former case aerogel layer
// with ref. index equal to the number is assumed, in the latter - material with name given by string.

class AddLayerCommand : public G4UIcommand {
public:
	AddLayerCommand(const char *theCommandPath,G4UImessenger *theMessenger);

	static void GetNewLayerSpecification(const char* paramString,G4String& layerSpec,G4double& t);
	//Parse string of command arguments and get material name/refractive index of aerogel as the
	// first argument and a thickness value specified in the last two arguments
};

#endif
