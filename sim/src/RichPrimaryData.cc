
#include <cmath>
#ifndef __CINT__
# warning "!!!__CINT__ not defined"
# include "G4ParticleDefinition.hh"
#endif
#include "RichPrimaryData.hh"

ClassImp(RichPrimaryData)

RichPrimaryData::~RichPrimaryData()
{}

RichPrimaryData::RichPrimaryData(const G4ParticleDefinition *part,Double_t p,Double_t theta,Double_t phi) :
	TObject(), fP(p), fTheta(theta), fPhi(phi)
{
	SetPrimaryData(part,p,theta,phi);
}

void RichPrimaryData::SetPrimaryData(const G4ParticleDefinition *part,Double_t p,Double_t theta,Double_t phi)
{
#ifndef __CINT__
	if( part==0 ) return;

    fPdgCode = part->GetPDGEncoding();
	fName = part->GetParticleName();
	fCharge = (int)part->GetPDGCharge();
	fMass = part->GetPDGMass()/GeV;
	fE = sqrt(fP*fP+fMass*fMass);
	fBeta = fP/fE;
#endif
}

