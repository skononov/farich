#include "RichAnalysisMessenger.hh"

#include "RichVAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

RichAnalysisMessenger::RichAnalysisMessenger(RichVAnalysisManager *
	                                                   analysisManager) :
	richAnalysis(analysisManager)
{

	RichAnalysisDir = new G4UIdirectory("/analysis/");
	RichAnalysisDir->SetGuidance("analysis control.");

	setInteractiveCommand = new G4UIcmdWithAString("/analysis/interactive", this);
	setInteractiveCommand->SetGuidance("Set interactive session on/off");
	setInteractiveCommand->SetParameterName("choice",true);
	setInteractiveCommand->SetDefaultValue("on");
	setInteractiveCommand->SetCandidates("on off");
	setInteractiveCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

	saveHitTreeCommand = new G4UIcmdWithABool("/analysis/saveHitTree",this);
	saveHitTreeCommand->SetGuidance("Set/Unset save hit tree flag.");
	saveHitTreeCommand->SetGuidance("If 'true' (default) the old hit tree will be written to root file.");
	saveHitTreeCommand->SetParameterName("choice",true);
	saveHitTreeCommand->SetDefaultValue(true);
	saveHitTreeCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

	saveEventTreeCommand = new G4UIcmdWithABool("/analysis/saveEventTree",this);
	saveEventTreeCommand->SetGuidance("Set/Unset save event tree flag.");
	saveEventTreeCommand->SetGuidance("If 'true' the new event tree will be written to root file.");
	saveEventTreeCommand->SetParameterName("choice",true);
	saveEventTreeCommand->SetDefaultValue(false);
	saveEventTreeCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

	runCintCommand = new G4UIcommand("/analysis/cint", this);
	runCintCommand->SetGuidance("Start ROOT C++ interpretator");
	runCintCommand->AvailableForStates(G4State_Idle);

	outputFileCommand = new G4UIcmdWithAString("/analysis/outputFile", this);
	outputFileCommand->SetGuidance("Specify the name of the output file.");
	outputFileCommand->SetGuidance("default: rich.root");
	outputFileCommand->SetParameterName("file.root",false);
	outputFileCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

   	macroFileCommand = new G4UIcmdWithAString("/analysis/macroFile", this);
	macroFileCommand->SetGuidance("Name of the macro file executed after each run.");
	macroFileCommand->SetGuidance("default: run.C");
	macroFileCommand->SetParameterName("macro_file",true);
	macroFileCommand->SetDefaultValue("run.C");
	macroFileCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

	betaDataFileCommand = new G4UIcmdWithAString("/analysis/betaDataFile", this);
	betaDataFileCommand->SetGuidance("Name of the file where beta vs radius data are stored.");
	betaDataFileCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

	verboseCommand = new G4UIcmdWithAnInteger("/analysis/verbose", this);
	verboseCommand->SetParameterName("verbosity",false);

	weighPixelsCommand = new G4UIcmdWithABool("/analysis/weighPixels",this);
	weighPixelsCommand->SetGuidance("Set/Unset weigh pixels flag.");
	weighPixelsCommand->SetGuidance("If 'true' weigh pixels by number of photoelectrons.");
	weighPixelsCommand->SetParameterName("choice",true);
	weighPixelsCommand->SetDefaultValue(false);
	weighPixelsCommand->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
}

RichAnalysisMessenger::~RichAnalysisMessenger()
{
	delete RichAnalysisDir;
}

void RichAnalysisMessenger::SetNewValue(G4UIcommand * command,
													G4String newValue)
{
	if (command==setInteractiveCommand) {
		G4bool b = (newValue=="on");
		richAnalysis->SetInteractive(b);
	}
	if (command==saveHitTreeCommand)
		richAnalysis->SetHitTreeStoring(saveHitTreeCommand->GetNewBoolValue(newValue));
	else if (command==saveEventTreeCommand)
		richAnalysis->SetEventTreeStoring(saveEventTreeCommand->GetNewBoolValue(newValue));
	else if (command==runCintCommand)
        richAnalysis->RunCint();
	else if (command==outputFileCommand)
		richAnalysis->SetOutputFileName(newValue);
	else if (command==macroFileCommand)
		richAnalysis->SetMacroFileName(newValue);
	else if (command==betaDataFileCommand)
		richAnalysis->SetBetaDataFileName(newValue);
	else if (command==verboseCommand)
		richAnalysis->SetVerboseLevel(verboseCommand->GetNewIntValue(newValue));
	else if (command==weighPixelsCommand)
		richAnalysis->SetWeighPixels(weighPixelsCommand->GetNewBoolValue(newValue));
}

/*
G4String RichAnalysisMessenger::GetCurrentValue(G4UIcommand * command)
{
	if (command==setInteractiveCommand)
		return G4String(richAnalysis->IsInteractive()?"on":"off");
	else if (command==outputFileCommand)
		return richAnalysis->GetOutputFileName();
	else if (command==macroFileCommand)
		return richAnalysis->GetMacroFileName();
	else if (command==verboseCommand)
		return verboseCommand->ConvertToString(richAnalysis->GetVerboseLevel());
	return G4String();
}
*/
