#ifndef RichFieldSetup_H
# define RichFieldSetup_H

# include "G4MagneticField.hh"
# include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class RichFieldMessenger;

class RichFieldSetup
{
  public:
	RichFieldSetup(G4ThreeVector);	//  The value of the field
	RichFieldSetup();			//  A zero field

	~RichFieldSetup();

	void SetStepperType(G4int i) { fStepperType = i; }

	void SetStepper();

	void SetMinStep(G4double s) { fMinStep = s; }

	void InitialiseAll();		//  Set parameters and call method below
	void CreateStepperAndChordFinder();

	void SetFieldValue(G4ThreeVector fieldVector) {
		((G4UniformMagField*)fMagneticField)->SetFieldValue(fieldVector);
	}

	G4ThreeVector GetFieldValue() const {
		return ((G4UniformMagField*)fMagneticField)->GetConstantFieldValue();
	}

  protected:

	// Find the global Field Manager

	G4FieldManager * GetGlobalFieldManager();	// static

	G4FieldManager *fFieldManager;
	G4ChordFinder *fChordFinder;
	G4Mag_UsualEqRhs *fEquation;
	G4MagneticField *fMagneticField;

	G4MagIntegratorStepper *fStepper;
	G4int fStepperType;

	G4double fMinStep;

	RichFieldMessenger *fFieldMessenger;

};

#endif
