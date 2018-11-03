#ifndef RichPrimaryData_h
# define RichPrimaryData_h 1

# include "TObject.h"
# include "TString.h"

class G4ParticleDefinition;

class RichPrimaryData : public TObject
{
private:
	Int_t    fPdgCode; //PDG particle code
	Int_t    fCharge;  //charge in units of electron charge
	TString  fName;    //name of the particle
	Double_t fMass;    //mass of the particle, GeV/c^2
	Double_t fP;       //absolute momentum of the particle, GeV/c
	Double_t fTheta;   //polar angle of the particle, deg
	Double_t fPhi;     //azimuthal angle of the particle, deg
	Double_t fE;       //energy of the particle, GeV
    Double_t fBeta;    //v/c of the particle

public:
	RichPrimaryData() {}
	RichPrimaryData(const G4ParticleDefinition *part,Double_t p,Double_t theta,Double_t phi);
	virtual ~RichPrimaryData();

	Int_t    PdgCode() const { return fPdgCode; }
	Int_t    Charge() const  { return fCharge; }
	TString  Name() const { return fName; }
	Double_t Mass() const { return fMass; }
	Double_t P() const { return fP; }
	Double_t E() const { return fE; }
	Double_t Beta() const { return fBeta; }
	Double_t Theta() const { return fTheta; }
	Double_t Phi() const { return fPhi; }

	void SetPrimaryData(const G4ParticleDefinition *part,Double_t p,Double_t theta,Double_t phi);

	ClassDefNV(RichPrimaryData,4)
};

#endif
