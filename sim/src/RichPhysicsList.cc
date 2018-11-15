
#include <iomanip>

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserPhysicsListMessenger.hh"
#include "G4UImanager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4UnitsTable.hh"

#include "G4Decay.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4MultipleScattering.hh"
#include "G4EmProcessOptions.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4SystemOfUnits.hh"

#include "RichCerenkov.hh"
#include "RichStepCut.hh"
#include "RichPhysicsList.hh"


RichPhysicsList::RichPhysicsList() :
	G4VUserPhysicsList()
{
	G4cout << "Now define the physics List\n";

	defaultCutValue = 0.1 * mm;
	constructed     = false;
	maxChargedStep  = 1 * mm;

	// pointer to the particle table
	theParticleTable = G4ParticleTable::GetParticleTable();
	theParticleIterator = theParticleTable->GetIterator();
}

RichPhysicsList::~RichPhysicsList()
{
}

void RichPhysicsList::ConstructParticle()
{
	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program.

	if (constructed) return;

	ConstructBosons();
	ConstructLeptons();
	ConstructMesons();
	ConstructBaryons();
}

void RichPhysicsList::ConstructBosons()
{
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
	G4ChargedGeantino::ChargedGeantinoDefinition();

	// gamma
	G4Gamma::GammaDefinition();

	// optical photon
	G4OpticalPhoton::OpticalPhotonDefinition();
}

void RichPhysicsList::ConstructLeptons()
{
	// leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void RichPhysicsList::ConstructMesons()
{
	//  mesons (charged)
	G4PionPlus::PionPlusDefinition();
	G4PionMinus::PionMinusDefinition();
	G4KaonPlus::KaonPlusDefinition();
	G4KaonMinus::KaonMinusDefinition();
}

void RichPhysicsList::ConstructBaryons()
{
//  barions
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
	G4Neutron::NeutronDefinition();
	G4AntiNeutron::AntiNeutronDefinition();
}

void RichPhysicsList::ConstructProcess()
{
	if (constructed) return;

	AddTransportation();
	ConstructGeneral();
	ConstructEM();
	ConstructOp();

	constructed = true;
}

void RichPhysicsList::ConstructGeneral()
{
	G4Decay *theDecayProcess = new G4Decay();
	theParticleIterator->reset();

#ifdef SECONDARIES_ON
	while ((*theParticleIterator) ()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		if (theDecayProcess->IsApplicable(*particle)) {
			pmanager->AddDiscreteProcess(theDecayProcess);
		}
	}
#endif
}

void RichPhysicsList::ConstructEM()
{
	theParticleIterator->reset();
	G4cout << " Now creating EM processes" << G4endl;

	RichStepCut* theStepCut = new RichStepCut;

#ifdef MS_ON
	G4cout << "  Multiple scattering process will be created" << G4endl;
#endif
#ifdef SECONDARIES_ON
	G4cout << "  Processes producing secandaries will be created" << G4endl;
#endif
	while ((*theParticleIterator) ()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

#ifdef MS_ON
		if ((particle->GetPDGCharge() != 0.0) &&
			(particle->GetParticleName() != "chargedgeantino")) {
			// all others charged particles except geantino
			pmanager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
		}
#endif

#ifdef SECONDARIES_ON
		if (particleName == "gamma") {
			// gamma
			// Construct processes for gamma
			pmanager->AddDiscreteProcess(new G4GammaConversion());
			pmanager->AddDiscreteProcess(new G4ComptonScattering());
			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

		} else if (particleName == "e-") {
			//electron
			//Construct processes for electron
//			pmanager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
			pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
			pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);

		} else if (particleName == "e+") {
			//positron
			// Construct processes for positron
//			pmanager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
			pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
			pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);
			pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);

		} else if (particleName == "mu+" || particleName == "mu-") {
			//muon
			// Construct processes for muon
//			pmanager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
			pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung(), -1, -1, 3);
			pmanager->AddProcess(new G4MuPairProduction(), -1, -1, 4);

		} else {
			if ((!particle->IsShortLived()) &&
				(particle->GetPDGCharge() != 0.0) &&
				(particle->GetParticleName() != "chargedgeantino")) {
				// all others charged particles except geantino
//				pmanager->AddProcess(new G4MultipleScattering(), -1, 1, 1);
				pmanager->AddProcess(new G4hIonisation(), -1, 2, 2);
			}
		}
#endif //SECONDARIES_ON

		if (theStepCut->IsApplicable(*particle) && !particle->IsShortLived()) {
			pmanager->AddDiscreteProcess(theStepCut);
			pmanager->SetProcessOrderingToLast(theStepCut,idxPostStep);
			theStepCut->SetMaxStep(maxChargedStep);
		}
	}
}

void RichPhysicsList::ConstructOp()
{
	G4cout << " Now creating Optical processes" << G4endl;
	G4OpAbsorption *theAbsorptionProcess = new G4OpAbsorption();
	G4OpRayleigh *theRayleighScatteringProcess = new G4OpRayleigh();
	G4OpBoundaryProcess *theBoundaryProcess = new G4OpBoundaryProcess();
	RichCerenkov *theCerenkovProcess = new RichCerenkov();

	theCerenkovProcess->SetVerboseLevel(0);
	theAbsorptionProcess->SetVerboseLevel(0);
	theRayleighScatteringProcess->SetVerboseLevel(0);
	theBoundaryProcess->SetVerboseLevel(0);

	G4int MaxNumPhotons = 300;

	theCerenkovProcess->SetTrackSecondariesFirst(true);
	theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

	G4OpticalSurfaceModel themodel = unified;
	theBoundaryProcess->SetModel(themodel);

	theParticleIterator->reset();
	while ((*theParticleIterator) ()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		if (theCerenkovProcess->IsApplicable(*particle)) {
			pmanager->AddContinuousProcess(theCerenkovProcess);
		}
		if (particleName == "opticalphoton") {
			pmanager->AddDiscreteProcess(theAbsorptionProcess);
			pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
			pmanager->AddDiscreteProcess(theBoundaryProcess);
		}
	}
}

void RichPhysicsList::SetCuts()
{
	if (verboseLevel > 1) {
		G4cout << "RichPhysicsList::SetCuts: ";
	}
	//  " G4VUserPhysicsList::SetCutsWithDefault" method sets
	//   the default cut value for all particle types
	SetCutsWithDefault();
}

