#include "RichParameters.hh"

using namespace std;

vector<G4double> defPhotMomVector = InitializePhotonMomentumVector();

G4double (*AerogelRefIndex)(G4double, G4double, G4bool) = AerogelRefIndexQuartz;

vector<G4double> InitializePhotonMomentumVector() {
	G4double PhotonEnergyStep=(kPhotonMaxEnergy-kPhotonMinEnergy)/kNPhotWLbins;
	vector<G4double> PhotMomVect(kNPhotWLbins+1);
	for (G4int ibin=0; ibin<=kNPhotWLbins; ibin++) {
		PhotMomVect[ibin]=kPhotonMinEnergy+PhotonEnergyStep*ibin;
	}
	return PhotMomVect;
}

G4double AerogelRefIndexConstant(G4double ri,G4double PhotMom,G4bool getRef) {
	return ri;
}

G4double AerogelRefIndexQuartz(G4double ri,G4double PhotMom,G4bool getRef) {
	G4double qzri=FusedSilicaRefIndex(PhotMom);
	if( getRef ) //ri is given at PhotMom, return ri at kReferencePhotMom
		return sqrt( 1 + (ri*ri-1)/(qzri*qzri-1)*(kReferenceQuartzRefInd*kReferenceQuartzRefInd-1) );
	else //ri is given at kReferencePhotMom, return ri at PhotMom
		return sqrt( 1 + (ri*ri-1)/(kReferenceQuartzRefInd*kReferenceQuartzRefInd-1)*(qzri*qzri-1) );
}

G4double AerogelRefIndexLHCb(G4double ri,G4double PhotMom,G4bool getRef) {
	G4double wl = kPhotMomWaveConv/(PhotMom/eV); //[nm]
	G4double ri2m1_lhcb = LHCb_a0/(1-LHCb_wl0sqr/(wl*wl));
	if( getRef ) //ri is given at PhotMom, return ri at kReferencePhotMom
		return sqrt( 1 + (ri*ri-1)/ri2m1_lhcb*LHCb_RI2m1ref );
	else //ri is given at kReferencePhotMom, return ri at PhotMom
		return sqrt( 1 + (ri*ri-1)/LHCb_RI2m1ref*ri2m1_lhcb );
}

