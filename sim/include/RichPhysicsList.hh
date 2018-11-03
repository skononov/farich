#ifndef RichPhysicsList_h
# define RichPhysicsList_h 1

# include "globals.hh"
# include "G4VUserPhysicsList.hh"
# include "G4ParticleTable

class RichPhysicsList : public G4VUserPhysicsList {
public:
	RichPhysicsList();
	virtual ~RichPhysicsList();

protected:
	// Construct particles and processes
	virtual void ConstructParticle();
	virtual void ConstructProcess();

	//
	virtual void SetCuts();

protected:
	// these methods Construct particles
	virtual void ConstructBosons();
	virtual void ConstructLeptons();
	virtual void ConstructMesons();
	virtual void ConstructBaryons();

protected:
	// these methods Construct physics processes and register them
	virtual void ConstructGeneral();
	virtual void ConstructEM();
	virtual void ConstructOp();

private:
	// the particle table has the complete List of existing particle types
	G4ParticleTable *theParticleTable;
	G4ParticleTable::G4PTblDicIterator *theParticleIterator;

	G4bool constructed;

	G4double maxChargedStep;
};

#endif
