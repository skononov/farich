#include "RichFieldSetup.hh"
#include "RichFieldMessenger.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

RichFieldSetup::RichFieldSetup(G4ThreeVector fieldVector)
:	fChordFinder(0), fStepper(0)
{
	fMagneticField = new G4UniformMagField(fieldVector);
	G4cout << " RichFieldSetup: magnetic field set to Uniform( "
		<< fieldVector << " ) " << G4endl;
	InitialiseAll();
}

RichFieldSetup::RichFieldSetup()
:	fChordFinder(0), fStepper(0)
{
	fMagneticField = new G4UniformMagField(G4ThreeVector(0,0,0));
	G4cout << " RichFieldSetup: magnetic field set to zero" << G4endl;
	InitialiseAll();
}

void RichFieldSetup::InitialiseAll()
{
	fFieldMessenger = new RichFieldMessenger(this);

	fEquation = new G4Mag_UsualEqRhs(fMagneticField);

	fMinStep = 0.1 * mm;		// minimal step of 0.1 mm is default

	fStepperType = 5;			// HelixExplicitEuler is default stepper

	fFieldManager = G4TransportationManager::GetTransportationManager()
		->GetFieldManager();

	CreateStepperAndChordFinder();
}

////////////////////////////////////////////////////////////////////////////////

RichFieldSetup::~RichFieldSetup()
{
	if(fMagneticField)
		delete fMagneticField;
	if(fChordFinder)
		delete fChordFinder;
	if(fStepper)
		delete fStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//

void RichFieldSetup::CreateStepperAndChordFinder()
{
	SetStepper();
	G4cout << "The minimal step is equal to " << fMinStep / mm << " mm" << G4endl;

	G4Field* field=fMagneticField;

	G4ThreeVector fieldVector = ((G4UniformMagField*)fMagneticField)->GetConstantFieldValue();

	G4cout << "Setting field vector to " << fieldVector / gauss << " gauss " << G4endl;

	if( fieldVector.mag() == 0. ) {
		// If the new field's value is Zero, signal it as below
		//   so that it is not used for propagation.
		field=0;
	}

	fFieldManager->SetDetectorField(field);
	fEquation->SetFieldObj(field);

	if(fChordFinder)
		delete fChordFinder;

	fChordFinder = new G4ChordFinder((G4MagneticField*)field, fMinStep, fStepper);

	fFieldManager->SetChordFinder(fChordFinder);

	return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void RichFieldSetup::SetStepper()
{
	if(fStepper)
		delete fStepper;

	switch (fStepperType) {
		case 0:
			fStepper = new G4ExplicitEuler(fEquation);
			G4cout << "G4ExplicitEuler is called" << G4endl;
			break;
		case 1:
			fStepper = new G4ImplicitEuler(fEquation);
			G4cout << "G4ImplicitEuler is called" << G4endl;
			break;
		case 2:
			fStepper = new G4SimpleRunge(fEquation);
			G4cout << "G4SimpleRunge is called" << G4endl;
			break;
		case 3:
			fStepper = new G4SimpleHeum(fEquation);
			G4cout << "G4SimpleHeum is called" << G4endl;
			break;
		case 4:
			fStepper = new G4ClassicalRK4(fEquation);
			G4cout << "G4ClassicalRK4 is called" << G4endl;
			break;
		case 5:
			fStepper = new G4HelixExplicitEuler(fEquation);
			G4cout << "G4HelixExplicitEuler (default) is called" << G4endl;
			break;
		case 6:
			fStepper = new G4HelixImplicitEuler(fEquation);
			G4cout << "G4HelixImplicitEuler is called" << G4endl;
			break;
		case 7:
			fStepper = new G4HelixSimpleRunge(fEquation);
			G4cout << "G4HelixSimpleRunge is called" << G4endl;
			break;
		case 8:
			fStepper = new G4CashKarpRKF45(fEquation);
			G4cout << "G4CashKarpRKF45 is called" << G4endl;
			break;
		case 9:
			fStepper = new G4RKG3_Stepper(fEquation);
			G4cout << "G4RKG3_Stepper is called" << G4endl;
			break;
		default:
			fStepper = 0;
	}
	return;
}

