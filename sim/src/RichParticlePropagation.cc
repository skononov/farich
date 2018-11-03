#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "G4SystemOfUnits.hh"

#include "RichParticlePropagation.hh"
#include "RichGlobals.hh"

const Double_t RichParticlePropagation::approx_limit = 0.005;
const Double_t RichParticlePropagation::delta = 1e-6;

RichParticlePropagation::RichParticlePropagation(const G4ThreeVector &field, const G4ThreeVector &startPosition, const G4ThreeVector &startMomentum, G4double charge) :
    inited(false), straight(true), funcZshift(NULL)
{
    Initialize(field,startPosition,startMomentum,charge);
}

void RichParticlePropagation::Initialize(const G4ThreeVector &field, const G4ThreeVector &startPosition, const G4ThreeVector &startMomentum, G4double charge)
{
	G4double Bmag = field.mag();
    position = startPosition;
	tangent = startMomentum.unit();
	if( tangent.z()<=0.0 ) {
		G4cerr << "RichParticlePropagation::Initialize(): Particle does not move forward in Z." << G4endl;
        return;
	}
	if( Bmag<=kTolerance*gauss )
		straight = true;
	else {
		straight = false;
        //magnetic field direction
		fieldDir = field.unit();
        //fractions of tangent vector that are parallel and perpendicular to field direction
		vPar = tangent.dot(fieldDir);
		vPerp = sqrt(1-vPar*vPar);
        //if motion is along the field extrapolate a straight track
		if( vPerp<kTolerance )
			straight = true;
		else {
	        //unit vector along tangent component perpendicular to field direction
			perpTangent = (tangent-vPar*fieldDir)/vPerp;
	        //curvature unit vector (force direction)
			curvature = charge*tangent.cross(fieldDir)/vPerp;
			//Lorentz radius
			radius = fabs(startMomentum.mag()/eV*vPerp/(300.*charge*Bmag/gauss))*cm;
            funcZshift = new TF1("funcZshift",RichParticlePropagation::EvalZshift,0.,1.,6);
			funcZshift->SetParameters(radius, perpTangent.z(), curvature.z(), fieldDir.z(), vPar, vPerp);
            hperiod = 2*TMath::Pi()*radius/vPerp;
			hmax = funcZshift->GetMaximumX(0.,hperiod);
		}
	}
    inited = true;
}

G4bool RichParticlePropagation::Propagate(G4double z0, G4ThreeVector& endPosition, G4ThreeVector& endDirection)
{
	if( !inited ) {
		G4cerr << "RichParticlePropagation::Propagate(): Initialize first!" << G4endl;
        return false;
	}
	if( z0<position.z() ) {
		G4cerr << "RichParticlePropagation::Propagate(): Particle moves in negative direction on Z." << G4endl;
		return false;
	}
	if( straight ) {
		endPosition = position +  tangent * (z0-position.z())/tangent.z();
		endDirection = tangent;
        return true;
	}
	G4double dz = z0-position.z(), H = hmax;
	if( funcZshift->Eval(H)<dz ) {
		G4double dzper = funcZshift->Eval(hperiod);
		if( dzper<kTolerance ) {
			G4cerr << "RichParticlePropagation::Propagate(): Particle can not reach z=" << z0 << G4endl;
			return false;
		}
        H = hperiod/dzper*dz;
	}
    G4double h = funcZshift->GetX(dz,0.,H);
	//check for h-value sanity
	if( h<kTolerance && fabs(funcZshift->Eval(h)-dz)>delta*mm ) {
		G4cerr << "RichParticlePropagation::Propagate(): Error in calculation of z value." << G4endl;
		return false;
	}
    //look for a smaller root
	if( funcZshift->Derivative(h)<0 ) {
		G4double h2 = funcZshift->GetX(dz,0.,h*0.999);
		if( h2<h*0.999 && fabs(funcZshift->Eval(h2)-dz)<=delta*mm )
			h = h2;
	}
	//evaluate end position and direction
	G4double Phi = h*vPerp/radius, SinPhi, CosPhi;
	if( fabs(Phi) > approx_limit ) {
		SinPhi = sin(Phi);
		CosPhi = cos(Phi);
	} else {
		G4double Phi2 = Phi*Phi;
		G4double Phi3 = Phi2 * Phi;
		G4double Phi4 = Phi2 * Phi2;
		SinPhi = Phi - 1.0/6.0 * Phi3;
		CosPhi = 1 - 0.5 * Phi2 + 1.0/24.0 * Phi4;
	}
	endPosition = position + radius*(SinPhi*perpTangent+(1-CosPhi)*curvature) + h*vPar*fieldDir;
	endDirection = vPerp*(perpTangent*CosPhi + curvature*SinPhi) + vPar*fieldDir;

	return true;
}

Double_t RichParticlePropagation::EvalZshift(Double_t *xx,Double_t *par)
{
	Double_t &h = *xx; //natural parameter of the trajectory (helix length)
	if( h==0. ) return 0;
	Double_t &R = par[0]; //Lorentz radius
	Double_t &zPtan = par[1], &zCurv = par[2], &zB = par[3]; //z-components of the perp. tangent, curvature and field unit vectors
    Double_t &vPar = par[4], &vPerp = par[5]; //parallel and perpendicular to B fraction of tangent

	Double_t Phi = h*vPerp/R, CosPhi, SinPhi;

    if( fabs(Phi) > approx_limit ) {
    	SinPhi = sin(Phi);
       	CosPhi = cos(Phi);
    } else {
		G4double Phi2 = Phi*Phi;
		G4double Phi3 = Phi2 * Phi;
		G4double Phi4 = Phi2 * Phi2;
		SinPhi = Phi - 1.0/6.0 * Phi3;
		CosPhi = 1 - 0.5 * Phi2 + 1.0/24.0 * Phi4;
    }

    return R*(SinPhi*zPtan + (1-CosPhi)*zCurv) + h*vPar*zB;
}

void RichParticlePropagation::Draw() const
{
	TApplication app("",0,0,0,0);

	TCanvas c("c1","c1",500,540);
	c.Draw();

	funcZshift->SetRange(0.,2*hperiod);
    funcZshift->Draw();

	app.Run();
}

