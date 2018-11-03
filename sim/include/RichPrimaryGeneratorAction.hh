#ifndef RichPrimaryGeneratorAction_h
# define RichPrimaryGeneratorAction_h 1

# include <cmath>
# include "G4VUserPrimaryGeneratorAction.hh"
# include "globals.hh"

# include "G4ParticleGun.hh"

class G4Event;
class G4ParticleDefinition;
class RichDetectorConstruction;
class RichPrimaryGeneratorMessenger;
class TF1;

class RichPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
	RichPrimaryGeneratorAction(RichDetectorConstruction*);
	~RichPrimaryGeneratorAction();

public:

	void GeneratePrimaries(G4Event *);	//concrete implementation
	//of base class virtual method

	void InitializeRun(); //set primary vertex position from detector geometry

	void Print() const;

	G4ParticleGun* GetGun() const { return particleGun; }

	const G4ParticleDefinition* GetParticle() const {
		return particleGun->GetParticleDefinition();
	}

	const G4String& GetParticleName() const {
		return particleGun->GetParticleDefinition()->GetParticleName();
	}

	G4double GetMass() const {
		return particleGun->GetParticleDefinition()->GetPDGMass();
	}

    //Charge in units of elementary charge
	G4double GetCharge() const {
		return particleGun->GetParticleDefinition()->GetPDGCharge();
	}

	G4double GetKineticEnergy() const { return particleGun->GetParticleEnergy(); }
	G4double GetEnergy() const { return GetKineticEnergy()+GetMass(); }

	G4double GetTotalMomentum() const {
		G4double E = GetEnergy();
		G4double M = GetMass();
		return sqrt(E*E-M*M);
	}

	G4double GetBeta() const {
		return GetTotalMomentum()/GetEnergy();
	}

	G4ParticleMomentum GetMomentumDirection() const {
		return particleGun->GetParticleMomentumDirection();
	}

	G4ThreeVector GetPosition() const {
		return particleGun->GetParticlePosition();
	}

	G4bool IsVertexFixed() const {
		return fixVertexPosition;
	}

	void SetPrimaryName(const G4String&);
	//set particle direction in spherical coordinates,
	void SetTheta(G4double);
	void SetPhi(G4double);
	void SetTotalMomentum(G4double);
	void FixVertexPosition(G4bool fix=true) { fixVertexPosition=fix; }

private:

	G4ParticleGun *particleGun;	//pointer a to G4 service class
	RichDetectorConstruction *richDetector;
	G4bool fixVertexPosition; //fix primary vertex position on radiator input face

	TF1 *fMomSpectrum; //spectrum of electrons for BGE_ON option

	RichPrimaryGeneratorMessenger *gunMessenger;	//messenger of this class
};

#endif
