#include "G4ThreeVector.hh"
#include "G4Timer.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Electron.hh"

#include "TF1.h"
#include "TGraph.h"

#include "RichPrimaryGeneratorAction.hh"
#include "RichParameters.hh"
#include "RichVAnalysisManager.hh"
#include "RichPrimaryGeneratorMessenger.hh"
#include "RichDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"


#ifdef BGE_ON
static TGraph* gmsp=new TGraph("/home/skononov/geant4/farich/superb/ebckgmom.dat");
static const G4double me=G4Electron::Electron()->GetPDGMass();

Double_t fmsp(Double_t *x,Double_t *par)
{
	if( gmsp->IsZombie() ) return 0.;

	Double_t p=*x;
	Double_t val=gmsp->Eval(p);

	if( p>3. )
		val=gmsp->Eval(3.)*std::exp(-1.9*p);

	return val<0.?0.:val;
}
#endif

RichPrimaryGeneratorAction::RichPrimaryGeneratorAction(
	RichDetectorConstruction* det) :
    G4VUserPrimaryGeneratorAction(),
	richDetector(det)
{
	G4int n_particle = 1;
	particleGun = new G4ParticleGun(n_particle);

	//create a messenger for this class
	gunMessenger = new RichPrimaryGeneratorMessenger(this);

	//default particle kinematic
	G4String thePrimaryParticleName = "pi-";
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

	G4ParticleDefinition *particle
		= particleTable->FindParticle(thePrimaryParticleName);

	particleGun->SetParticleDefinition(particle);

	particleGun->SetParticleMomentum(G4ThreeVector(0.,0.,4.*GeV));

	InitializeRun();

	fixVertexPosition=false;

#ifdef BGE_ON
	particleGun->SetParticleDefinition(G4Electron::Electron());
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	particleGun->SetParticleEnergy(1.*GeV);
	fMomSpectrum=new TF1("fmsp",fmsp,0.,3.,0);
	fMomSpectrum->SetNpx(120);
#endif
}

RichPrimaryGeneratorAction::~RichPrimaryGeneratorAction()
{
#ifdef BGE_ON
	delete fMomSpectrum;
#endif
	delete particleGun;
	delete gunMessenger;
}

void RichPrimaryGeneratorAction::InitializeRun()
{
	particleGun->SetParticlePosition(G4ThreeVector(0,0,richDetector->GetRadiatorInSurfZ()-1*mm));
}

void RichPrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
	G4ThreeVector position(0,0,richDetector->GetRadiatorInSurfZ()-1*mm);

	G4double period=richDetector->GetPixelSpacing();

	//set primary position randomly over the square of size equal to pixel spacing
	if( !fixVertexPosition && period>0. ) {
		position.setX((G4UniformRand()-0.5)*period);
		position.setY((G4UniformRand()-0.5)*period);
	}

#ifdef BGE_ON
	G4double costh=G4UniformRand();
	G4double sinth=std::sqrt(1-costh*costh);
	G4double phi=2*pi*G4UniformRand();
	G4ThreeVector dir(sinth*std::cos(phi),sinth*std::sin(phi),costh);
	particleGun->SetParticleMomentumDirection(dir);
//	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	position.setZ(richDetector->GetRadiatorInSurfZ()+0.01*mm);
	G4double p=fMomSpectrum->GetRandom()*MeV;
	G4double kinetic_energy=std::sqrt(p*p+me*me)-me;
	particleGun->SetParticleEnergy(kinetic_energy);
#endif

	particleGun->SetParticlePosition(position);

	//this function is called at the beginning of event
	particleGun->GeneratePrimaryVertex(anEvent);
}

///////////////////////////////////////////////////////////////////////
//
//

void RichPrimaryGeneratorAction::SetPrimaryName(const G4String& particleName)
{
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	if (!particle) {
		G4cerr << " Particle " << particleName << " is unkown.\n";
		return;
	}
	G4double momentum=GetTotalMomentum();
	particleGun->SetParticleDefinition(particle);
	SetTotalMomentum(momentum);
}
void RichPrimaryGeneratorAction::SetTheta(G4double theta)
{
	G4ThreeVector direction = particleGun->GetParticleMomentumDirection();

	direction.setTheta(theta / radian);

	particleGun->SetParticleMomentumDirection(direction);

//	G4cout << " Direction of primary: Theta= "
//		   << direction.theta() / radian * degree << "\n";
}
void RichPrimaryGeneratorAction::SetPhi(G4double phi)
{
	G4ThreeVector direction = particleGun->GetParticleMomentumDirection();

	direction.setPhi(phi / radian);

	particleGun->SetParticleMomentumDirection(direction);

//	G4cout << " Direction of primary: Phi="
//		   << direction.phi() / radian * degree << ".\n";
}
void RichPrimaryGeneratorAction::SetTotalMomentum(G4double aMomentum)
{
	G4ThreeVector momentum = aMomentum*particleGun->GetParticleMomentumDirection();
	particleGun->SetParticleMomentum(momentum);
//	G4cout << " Total momentum of primary "
//		   << momentum.mag() / GeV << " GeV/c.\n";
}
void RichPrimaryGeneratorAction::Print() const
{
	G4ThreeVector direction=particleGun->GetParticleMomentumDirection();

	G4cout.precision(5);
	G4cout<<"Particle: "<<particleGun->GetParticleDefinition()->GetParticleName()<<"\n"
		  <<"Momentum="<<GetTotalMomentum()/GeV<<" GeV/c\n"
		  <<"Beta="<<GetBeta()<<"\n"
		  <<"Position: "<<(fixVertexPosition?"fixed":"random")<<"\n"
		  <<"Direction: theta="<<direction.theta()/degree
		  <<"deg, phi="<<direction.phi()/degree<<"deg"<<G4endl;
}

