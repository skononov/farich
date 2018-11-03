#include "RichFieldMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "RichGlobals.hh"
#include "RichFieldSetup.hh"

//////////////////////////////////////////////////////////////////////////////

RichFieldMessenger::RichFieldMessenger(RichFieldSetup * pEMfieldSetup)
:	fEMfieldSetup(pEMfieldSetup)
{
	FieldDir = new G4UIdirectory("/field/");
	FieldDir->SetGuidance("Rich field tracking control.");

	StepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType", this);
	StepperCmd->SetGuidance("Select stepper type for magnetic field");
	StepperCmd->SetParameterName("choice", true);
	StepperCmd->SetDefaultValue(4);
	StepperCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	UpdateCmd = new G4UIcmdWithoutParameter("/field/update", this);
	UpdateCmd->SetGuidance("Update field value and min step.");
	UpdateCmd->AvailableForStates(G4State_Idle);

	MagFieldCmd = new RichMagFieldCommand("/field/setField", this);

	MinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep", this);
	MinStepCmd->SetGuidance("Define minimal step");
	MinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
	MinStepCmd->SetParameterName("min_step", false, false);
	MinStepCmd->SetDefaultUnit("mm");
	MinStepCmd->AvailableForStates(G4State_Idle);
}

///////////////////////////////////////////////////////////////////////////////

RichFieldMessenger::~RichFieldMessenger()
{
	delete StepperCmd;
	delete MagFieldCmd;
	delete MinStepCmd;
	delete UpdateCmd;
	delete FieldDir;
}

////////////////////////////////////////////////////////////////////////////

void RichFieldMessenger::SetNewValue(G4UIcommand * command, G4String newValue)
{
	if(command == StepperCmd) {
		fEMfieldSetup->SetStepperType(StepperCmd->GetNewIntValue(newValue));
	}
	if(command == UpdateCmd) {
		fEMfieldSetup->CreateStepperAndChordFinder();
	}
	if(command == MagFieldCmd) {
		fEMfieldSetup->SetFieldValue(MagFieldCmd->GetFieldVector(newValue));
	}
	if(command == MinStepCmd) {
		fEMfieldSetup->SetMinStep(MinStepCmd->GetNewDoubleValue(newValue));
	}
}

/////////////////////////////////////////////////////////////////////////
//
//
#include <G4UIparameter.hh>
#include <G4Tokenizer.hh>

RichMagFieldCommand::RichMagFieldCommand(const char *theCommandPath,G4UImessenger *theMessenger) :
	G4UIcommand(theCommandPath,theMessenger)
{
	SetGuidance("/field/setField field [field_unit] [theta] [phi] [angle_unit]");
	SetGuidance("Define magnetic field. Magnetic field direction is given by angles.");

	G4UIparameter *fieldParameter=new G4UIparameter("field",'d',false);
	SetParameter(fieldParameter);

	G4UIparameter *fieldUnitParameter=new G4UIparameter("field_unit",'s',true);
	SetParameter(fieldUnitParameter);
	fieldUnitParameter->SetDefaultValue("T");
	fieldUnitParameter->SetParameterCandidates(UnitsList("Magnetic flux density"));

	G4UIparameter *thetaAngleParameter=new G4UIparameter("theta",'d',true);
	SetParameter(thetaAngleParameter);
	thetaAngleParameter->SetDefaultValue(0.0);

	G4UIparameter *phiAngleParameter=new G4UIparameter("phi",'d',true);
	SetParameter(phiAngleParameter);
	phiAngleParameter->SetDefaultValue(0.0);

	G4UIparameter *angleUnitParameter=new G4UIparameter("angle_unit",'s',true);
	SetParameter(angleUnitParameter);
	angleUnitParameter->SetDefaultValue("deg");
	angleUnitParameter->SetParameterCandidates(UnitsList("Angle"));
}

G4ThreeVector RichMagFieldCommand::GetFieldVector(const char* paramString)
{
	G4Tokenizer argTokenizer(paramString);

	G4String arg1=argTokenizer(), arg2=argTokenizer(),
		     arg3=argTokenizer(), arg4=argTokenizer(), arg5=argTokenizer();

	G4double fieldStrength = ConvertToDouble(arg1) * ValueOf(arg2);

	if( fieldStrength<kTolerance )
		return G4ThreeVector(0.,0.,0.);

	G4double angleUnit=ValueOf(arg5);
	G4double theta=ConvertToDouble(arg3) * angleUnit;
	G4double phi=ConvertToDouble(arg4) * angleUnit;

	G4double Bx=fieldStrength*sin(theta/rad)*cos(phi/rad),
		     By=fieldStrength*sin(theta/rad)*sin(phi/rad),
		     Bz=fieldStrength*cos(theta/rad);

	if( fabs(Bx)<kTolerance ) Bx = 0.;
	if( fabs(By)<kTolerance ) By = 0.;
	if( fabs(Bz)<kTolerance ) Bz = 0.;

	return G4ThreeVector(Bx,By,Bz);
}
