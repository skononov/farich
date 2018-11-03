#include "RichDetectorMessenger.hh"
#include "RichDetectorConstruction.hh"
#include "RichParameters.hh"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4Tokenizer.hh>

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichDetectorMessenger::RichDetectorMessenger(RichDetectorConstruction* detector) :
	richDetector(detector)
{
	char condition[100];

// Global setup
	richDir = new G4UIdirectory("/Rich/");
	richDir->SetGuidance("RICH detector UI command root.");

	verboseCmd = new G4UIcmdWithAnInteger("/Rich/verbose",this);
	verboseCmd->SetGuidance("Set detector construction verbosity level");
	verboseCmd->SetParameterName("level",false);

	defaultsCmd = new G4UIcmdWithoutParameter("/Rich/defaults",this);
	defaultsCmd->SetGuidance("Set all detector configuration values to defaults.");
	defaultsCmd->SetGuidance("(Update still required)");

	updateCmd = new G4UIcmdWithoutParameter("/Rich/update",this);
	updateCmd->SetGuidance("Run detector update.");

    printCmd = new G4UIcmdWithoutParameter("/Rich/print",this);
	printCmd->SetGuidance("Print out detector parameters.");

	overallDimensionCmd = new G4UIcmdWithADoubleAndUnit("/Rich/dimension", this);
	overallDimensionCmd->SetGuidance("Set overall dimension of aerogel and gap. Fixed-thickness mode only.");
	overallDimensionCmd->SetParameterName("dist",false);
	overallDimensionCmd->SetDefaultUnit("mm");
	overallDimensionCmd->SetRange("dist>0");

	proximityDistanceCmd = new G4UIcmdWithADoubleAndUnit("/Rich/proximity", this);
	proximityDistanceCmd->SetGuidance("Set distance between radiator and PMT array");
	proximityDistanceCmd->SetParameterName("dist",false);
	proximityDistanceCmd->SetDefaultUnit("mm");
	proximityDistanceCmd->SetRange("dist>0");

	selectModeCmd = new G4UIcmdWithAString("/Rich/mode",this);
	selectModeCmd->SetGuidance("Multi-layer radiator description mode.");
	selectModeCmd->SetGuidance("manual - manual radiator construction with addLayer command");
	selectModeCmd->SetGuidance("number - fix number of layers for single ring option");
	selectModeCmd->SetGuidance("thick  - fix total thickness for single ring option");
	selectModeCmd->SetGuidance("fixed - fix number of layers & total thickness for single ring option");
	selectModeCmd->SetGuidance("multiring - multi-ring radiator made of two series of aerogels");
	selectModeCmd->SetParameterName("choice",true);
    selectModeCmd->SetCandidates("manual number thick fixed multiring");

// Radiator setup
	radiatorDir = new G4UIdirectory("/Rich/radiator/");
	radiatorDir->SetGuidance("Radiator setup.");

	addLayerCmd = new AddLayerCommand("/Rich/radiator/addLayer",this);

	resetCmd = new G4UIcmdWithoutParameter("/Rich/radiator/reset",this);
	resetCmd->SetGuidance("Reset prefiously defined radiator layers.");

	normThickness1Cmd = new G4UIcmdWithADouble("/Rich/radiator/normThickness1", this);
    normThickness1Cmd->SetGuidance("Layer 1 normalized thickness. Fixed-thickness mode only.");
	normThickness1Cmd->SetParameterName("value",false);
	normThickness1Cmd->SetRange("value>0");

	layer1RefIndCmd = new G4UIcmdWithADouble("/Rich/radiator/index1", this);
	layer1RefIndCmd->SetGuidance("Set nominal refractive index of the first layer of radiator");
	layer1RefIndCmd->SetParameterName("n",false);
	layer1RefIndCmd->SetRange("n>=1.0");

	layer1ThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Rich/radiator/thickness1", this);
	layer1ThicknessCmd->SetGuidance("Set thickness of the first layer of radiator");
	layer1ThicknessCmd->SetParameterName("t",false);
	layer1ThicknessCmd->SetDefaultUnit("mm");
	layer1ThicknessCmd->SetRange("t>0");

	nLayersCmd = new G4UIcmdWithAnInteger("/Rich/radiator/nlayers", this);
	nLayersCmd->SetGuidance("Set number of radiator layers in radiator for fixed-number-of-layers mode");
	nLayersCmd->SetParameterName("N",false);
	sprintf(condition,"N>0 && N<=%d",kMaxNumberOfLayers);
	nLayersCmd->SetRange(condition);

	totalThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Rich/radiator/totalThickness", this);
	totalThicknessCmd->SetGuidance("Set thickness of the radiator");
	totalThicknessCmd->SetParameterName("t",false);
	totalThicknessCmd->SetDefaultUnit("mm");
	totalThicknessCmd->SetRange("t>0");

	betaOptimizedCmd = new G4UIcmdWithADouble("/Rich/radiator/betaOp",this);
	betaOptimizedCmd->SetGuidance("Set velocity value which radiator is to be optimized for");
	betaOptimizedCmd->SetParameterName("beta",false);
	betaOptimizedCmd->SetRange("beta>0 && beta<=1.0");

	series2RefIndCmd = new G4UIcmdWithADouble("/Rich/radiator/index2", this);
	series2RefIndCmd->SetGuidance("Set refractive index of a second series radiator layers. Multi-ring mode only.");
	series2RefIndCmd->SetParameterName("n",false);
	series2RefIndCmd->SetRange("n>=1.0");

	// Distortion setup
	distortThicknessCmd = new G4UIcmdWithADouble("/Rich/radiator/distortThickness",this);
	distortThicknessCmd->SetGuidance("Distort thickness of every layer by some relative magnitude");
	distortThicknessCmd->SetParameterName("eps",false);

	distortIndexCmd = new G4UIcmdWithADouble("/Rich/radiator/distortIndex",this);
	distortIndexCmd->SetGuidance("Distort ref. index of every layer by some relative density magnitude");
	distortIndexCmd->SetParameterName("eps",false);

	distortUniformityCmd = new G4UIcmdWithADouble("/Rich/radiator/distortUniformity",this);
	distortUniformityCmd->SetGuidance("Distort ref. index uniformity along beamline by some relative density magnitude");
	distortUniformityCmd->SetParameterName("eps",false);

	evenDistortionCmd = new G4UIcmdWithABool("/Rich/radiator/evenDistortion",this);
	evenDistortionCmd->SetGuidance("Set distortion to be even with layer number");
	evenDistortionCmd->SetParameterName("on",true);
	evenDistortionCmd->SetDefaultValue(true);

// Aerogel material optical properties
	scatteringLengthCmd = new G4UIcmdWithADoubleAndUnit("/Rich/radiator/scatterLength", this);
	scatteringLengthCmd->SetGuidance("Set Rayleigh scattering length in radiator at 400nm");
	scatteringLengthCmd->SetParameterName("value",true);
	scatteringLengthCmd->SetDefaultUnit("mm");
	scatteringLengthCmd->SetDefaultValue(defAerogelScatLength);
	scatteringLengthCmd->SetRange("value>0");

	absLengthDataFileCmd = new G4UIcmdWithAString("/Rich/radiator/absDataFile", this);
	absLengthDataFileCmd->SetGuidance("Set data file name for absorption length in radiator");
	absLengthDataFileCmd->SetParameterName("filename",true);
	absLengthDataFileCmd->SetDefaultValue(defAerogelAbsLenDataFile);

	boundaryLossCmd = new G4UIcmdWithADouble("/Rich/radiator/boundaryLoss", this);
    boundaryLossCmd->SetGuidance("Set loss ratio at radiator boundary");
	boundaryLossCmd->SetParameterName("value",true);
	boundaryLossCmd->SetDefaultValue(defAerogelSurfTrans);
	boundaryLossCmd->SetRange("value>=0 && value<=1");

	chromaticityCmd = new G4UIcmdWithAString("/Rich/radiator/chromaticity", this);
	chromaticityCmd->SetGuidance("Choice of dispersion model that determines dependence of refractive index on wavelength.");
	chromaticityCmd->SetGuidance("off        - no dispersion, constant refractive index;");
	chromaticityCmd->SetGuidance("on, quartz - default dispersion model as for quartz scaled by density to aerogel;");
	chromaticityCmd->SetGuidance("lhcb       - LHCb Sellmeier fit for n=1.03 scaled by density;");
	chromaticityCmd->SetParameterName("choice",true);
	chromaticityCmd->SetCandidates("off on quartz lhcb");
	chromaticityCmd->SetDefaultValue("quartz");

// Photodetector setup
	pmtDir = new G4UIdirectory("/Rich/pmt/");
	pmtDir->SetGuidance("PMT detector parameters.");

	qeDataFileCmd = new G4UIcmdWithAString("/Rich/pmt/qeDataFile", this);
	qeDataFileCmd->SetGuidance("Set data file name for spectral quantum efficiency");
	qeDataFileCmd->SetParameterName("name",true);
	qeDataFileCmd->SetDefaultValue(defPMTQEdataFile);

	detectionEffCmd = new G4UIcmdWithADouble("/Rich/pmt/detection", this);
	detectionEffCmd->SetGuidance("Set overall detection efficiency to PMT");
	detectionEffCmd->SetParameterName("efficiency",true);
	detectionEffCmd->SetDefaultValue(defDetectionEfficiency);
	detectionEffCmd->SetRange("efficiency>=0 && efficiency<=1");

	geomEffCmd = new G4UIcmdWithADouble("/Rich/pmt/geomEfficiency", this);
	geomEffCmd->SetGuidance("Set geometrical efficiency (active/total area fraction)");
	geomEffCmd->SetParameterName("efficiency",true);
	geomEffCmd->SetDefaultValue(defGeomEfficiency);
	geomEffCmd->SetRange("efficiency>=0 && efficiency<=1");

	pixelSizeCmd = new G4UIcmdWithADoubleAndUnit("/Rich/pmt/pixelSize", this);
	pixelSizeCmd->SetGuidance("Set pixel size of photodetector");
	pixelSizeCmd->SetParameterName("size",false);
	pixelSizeCmd->SetDefaultUnit("mm");
	pixelSizeCmd->SetRange("size>=0");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichDetectorMessenger::~RichDetectorMessenger()
{
	delete richDir;
	delete radiatorDir;
	delete pmtDir;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == verboseCmd)
		richDetector->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));
	else if (command == defaultsCmd)
        richDetector->SetDefaults();
	else if (command == updateCmd)
		richDetector->Update();
	else if (command == printCmd)
		richDetector->Print();
	else if (command == overallDimensionCmd)
		richDetector->SetOverallDimension(
			overallDimensionCmd->GetNewDoubleValue(newValue));
	else if (command == selectModeCmd) {
		if (newValue=="manual")
			richDetector->SetMLRmode(RichDetectorConstruction::manual);
		else if (newValue=="number")
			richDetector->SetMLRmode(RichDetectorConstruction::fixedNumberOfLayers);
		else if (newValue=="thick")
			richDetector->SetMLRmode(RichDetectorConstruction::fixedThickness);
		else if (newValue=="fixed")
			richDetector->SetMLRmode(RichDetectorConstruction::fixedNT);
		else if (newValue=="multiring")
			richDetector->SetMLRmode(RichDetectorConstruction::multiRing);
	}
	else if (command == addLayerCmd) {
		G4String spec;
		G4double t;
		addLayerCmd->GetNewLayerSpecification(newValue,spec,t);
		richDetector->AddLayer(spec,t);
	}
	else if (command == resetCmd)
		richDetector->ResetRadiator();
	else if (command == normThickness1Cmd)
		richDetector->SetThickness1Normalized(
			normThickness1Cmd->GetNewDoubleValue(newValue));
	else if (command == nLayersCmd)
		richDetector->SetNumberOfLayers(
			nLayersCmd->GetNewIntValue(newValue));
	else if (command == proximityDistanceCmd)
		richDetector->SetProximityDistance(
			proximityDistanceCmd->GetNewDoubleValue(newValue));
	else if (command == layer1RefIndCmd)
		richDetector->SetLayer1NomRefIndex(
			layer1RefIndCmd->GetNewDoubleValue(newValue));
	else if (command == layer1ThicknessCmd)
		richDetector->SetLayer1Thickness(
			layer1ThicknessCmd->GetNewDoubleValue(newValue));
	else if (command == totalThicknessCmd)
		richDetector->SetTotalThickness(
			totalThicknessCmd->GetNewDoubleValue(newValue));
	else if (command == betaOptimizedCmd)
		richDetector->SetBetaOptimized(
			betaOptimizedCmd->GetNewDoubleValue(newValue));
	else if (command == series2RefIndCmd)
		richDetector->SetSeries2NomRefIndex(
			series2RefIndCmd->GetNewDoubleValue(newValue));
	else if (command == distortThicknessCmd)
		richDetector->SetThicknessDistortion(
			distortThicknessCmd->GetNewDoubleValue(newValue));
	else if (command == distortIndexCmd)
		richDetector->SetDensityDistortion(
			distortIndexCmd->GetNewDoubleValue(newValue));
	else if (command == distortUniformityCmd)
		richDetector->SetUniformityDistortion(
			distortUniformityCmd->GetNewDoubleValue(newValue));
	else if (command == evenDistortionCmd)
		richDetector->SetEvenDistortion(
			evenDistortionCmd->GetNewBoolValue(newValue));
	else if (command == scatteringLengthCmd)
		richDetector->SetAerogelScatLength(
			scatteringLengthCmd->GetNewDoubleValue(newValue));
	else if (command == absLengthDataFileCmd) {
        if (newValue.contains('/'))
			richDetector->SetAerogelAbsorptionFile(newValue);
		else
			richDetector->SetAerogelAbsorptionFile(kInputDataDir+newValue);
	}
	else if (command == boundaryLossCmd)
		richDetector->SetAerogelBoundaryLoss(boundaryLossCmd->GetNewDoubleValue(newValue));
	else if (command == chromaticityCmd) {
		if (newValue=="off")
			richDetector->SetChromaticity(kNoDispersion);
		else if (newValue=="on" || newValue=="quartz")
			richDetector->SetChromaticity(kQuartzModel);
		else if (newValue=="lhcb")
			richDetector->SetChromaticity(kLHCbModel);
	}
	else if (command == qeDataFileCmd) {
        if (newValue.contains('/'))
			richDetector->SetQEfile(newValue);
		else
			richDetector->SetQEfile(kInputDataDir+newValue);
	}
	else if (command == detectionEffCmd)
		richDetector->SetDetectionEfficiency(
			detectionEffCmd->GetNewDoubleValue(newValue));
	else if (command == geomEffCmd)
		richDetector->SetGeomEfficiency(
			geomEffCmd->GetNewDoubleValue(newValue));
	else if (command==pixelSizeCmd)
		richDetector->SetPixelSize(pixelSizeCmd->GetNewDoubleValue(newValue));
}

AddLayerCommand::AddLayerCommand(const char *theCommandPath,G4UImessenger *theMessenger) :
	G4UIcommand(theCommandPath,theMessenger)
{
	SetGuidance("Add a radiator layer manually giving material name or aerogel index of refraction");
	SetGuidance("as the first argument, thickness and its unit as the second and third.");

	G4UIparameter *specParameter=new G4UIparameter("n",'s',false);
	SetParameter(specParameter);

	G4UIparameter *thicknessParameter=new G4UIparameter("t",'d',false);
	SetParameter(thicknessParameter);
	G4UIparameter *thicknessUnit=new G4UIparameter("unit",'s',true);
	SetParameter(thicknessUnit);
	thicknessUnit->SetDefaultValue("mm");
	thicknessUnit->SetParameterCandidates(UnitsList("Length"));
}

void AddLayerCommand::GetNewLayerSpecification(const char* paramString,G4String& layerSpec,G4double& t)
{
	G4Tokenizer argTokenizer(paramString);

	layerSpec = argTokenizer(); //first argument is a material name or index value of aerogel
	t = ConvertToDouble(G4String(argTokenizer())) * ValueOf(G4String(argTokenizer()));
}

