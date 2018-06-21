#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "Math/SpecFunc.h"

using namespace std;

const Double_t eulergamma=0.5772156649;
const Double_t Const = 2*TMath::Pi()/137.036*1e6; //Cherenkov emission coefficient assuming wavelength in nm and path in mm
const Double_t eps = 1e-10;

TGraph *gQE = 0;

Int_t Nlambda = 50, Nbins = 200; 
Double_t Lsc = 50; //nm

struct Result {
    Bool_t valid;
    Double_t npe, R, Ang, sigR, sigR_t, sigAng, sigAng_t;
    TH1* hrad;
};

/*Refractive index as function of wavelength
    for aerogel of the nominal ref. ind. at 400nm.
  Scaled from LHCb data:
  T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
  Eur. Phys. J. C 52 (2007) 759-764
*/
Double_t aerindex(Double_t n400,Double_t wl)
{
    const double LHCb_a0=0.05639, LHCb_wl0sqr = 83.22*83.22;
    double LHCb_RI2m1ref = LHCb_a0/(1-LHCb_wl0sqr/(400*400)); //(n**2-1) of LHCb aerogel at 400nm
    double ri2m1_lhcb = LHCb_a0/(1-LHCb_wl0sqr/(wl*wl)); //(n**2-1) of LHCb aerogel at wl

    return sqrt( 1 + (n400*n400-1)/LHCb_RI2m1ref*ri2m1_lhcb );
}

// Tangent of angle at point x for the photon emitted at point x0 for given radiator and wavelength
Double_t func_tan(Double_t *xx, Double_t *par)
{
    Double_t x = *xx;
    Double_t n_i = par[0], k = par[1], b2 = par[2], x0 = par[3], wl = par[4];
    Double_t n_wl = aerindex(n_i+k*x,wl), n0_wl = aerindex(n_i+k*x0,wl);
    
    return 1/sqrt(n_wl*n_wl/(n0_wl*n0_wl-1/b2)-1);
}

// Inverse cosine of angle at point x for the photon emitted at point x0 for given radiator and wavelength
Double_t func_path(Double_t *xx, Double_t *par)
{
    Double_t x = *xx;
    Double_t n_i = par[0], k = par[1], b2 = par[2], x0 = par[3], wl = par[4];
    Double_t n_wl = aerindex(n_i+k*x,wl), n0_wl = aerindex(n_i+k*x0,wl);
    
    return sqrt(1+1/(n_wl*n_wl/(n0_wl*n0_wl-1/b2)-1));
}

TF1* ftan = 0;
TF1* fpath = 0;

Double_t func_r(Double_t *xx, Double_t *par)
{
	Double_t x0 = xx[0];
	Double_t n_i = par[0], k = par[1], La = par[2], Lv = par[3], b2 = par[4], wl = par[5];
	Double_t n0_wl = aerindex(n_i+k*x0,wl);
	
	ftan->SetParameters(n_i,k,b2,x0,wl);

	return ftan->Integral(x0,La) + Lv/sqrt(1/(n0_wl*n0_wl-1/b2)-1);
}

Double_t func_i(Double_t *xx, Double_t *par)
{
	Double_t x0 = xx[0], wl = xx[1];
	Double_t n_i = par[0], k = par[1], La = par[2], b2 = par[3];
	Double_t n0_wl = aerindex(n_i+k*x0,wl);
	
	fpath->SetParameters(n_i,k,b2,x0,wl);
	Double_t path = fpath->Integral(x0,La);
	
	return Const * (1-1/b2/(n0_wl*n0_wl)) / (wl*wl) * exp(-path/Lsc*pow(400./wl,4));
}

Result farich_dIdR(const char* qefn, Double_t n_i, Double_t n_f, Double_t La, Double_t Ltot, Double_t b, Bool_t draw=kTRUE)
{
    Result res;
    res.valid = kFALSE;

    gQE = new TGraph(qefn);
    if( gQE->GetN()==0 ) return res;

	Double_t b2 = b*b, Lv = Ltot-La;
	Double_t k = (n_f-n_i)/La;
	Double_t wl_min = gQE->GetX()[0];
	Double_t wl_max = gQE->GetX()[gQE->GetN()-1];
	Double_t dwl = (wl_max-wl_min)/Nlambda;

    ftan = new TF1("ftan",func_tan,0,La,5);
    fpath = new TF1("fpath",func_path,0,La,5);

	TF1* frad = new TF1("frad",func_r,0.,La,6);
	frad->SetParNames("n_i","k","La","Lv","b2","wl");
	frad->SetParameters(n_i,k,La,Lv,b2,wl_min);

	TF2* fint = new TF2("fint",func_i,0.,La,wl_min,wl_max,4);
	fint->SetParNames("n_i","k","La","b2");
	fint->SetParameters(n_i,k,La,b2);

	frad->SetParameter(5,wl_max);
	Double_t Rmin0 = TMath::Max(0.,frad->GetMinimum()-2);
	frad->SetParameter(5,wl_min);
	Double_t Rmax0 = frad->GetMaximum()+2;

	Double_t *I = new Double_t[Nbins];
	for(int i=0; i<Nbins; i++) I[i] = 0.;
	Double_t dR = (Rmax0-Rmin0)/Nbins;
	
	for(int j=0; j<Nlambda; j++){
		Double_t wl = wl_min+dwl*(j+0.5);

		Double_t qe = gQE->Eval(wl)/100;
		frad->SetParameter(5,wl);
		Double_t x_ext = frad->GetMaximumX();
		Double_t Rmin = frad->GetMinimum(), Rmax = frad->GetMaximum();
		Bool_t monotonic = kFALSE;
		if( x_ext<eps || La-x_ext<eps )
			monotonic = kTRUE;

		for(int i=0; i<Nbins; i++){
			Double_t R0 = Rmin0+i*dR;
			Double_t R1 = R0+dR;

			if( (R0 < Rmin && R1 < Rmin) || (R0 > Rmax && R1 > Rmax) )
				continue;

            R0 = R0<Rmin?Rmin:R0;
            R1 = R1>Rmax?Rmax:R1;

			Double_t x1[2]={0,0}, x2[2]={0,0};
			if( monotonic ) {			
				x1[0] = frad->GetX(R0);
				x1[1] = frad->GetX(R1);
			} else {
				x1[0] = frad->GetX(R0,0.,x_ext);
				x1[1] = frad->GetX(R1,0.,x_ext);
				x2[0] = frad->GetX(R0,x_ext,La);
				x2[1] = frad->GetX(R1,x_ext,La);
			}

		    Double_t x0 = (x1[0]+x1[1])/2;
			I[i] += fabs(x1[1]-x1[0])/dR*(*fint)(x0,wl)*qe*dwl;
			if( !monotonic ) {
			    x0 = (x2[0]+x2[1])/2;
				I[i] += fabs(x2[1]-x2[0])/dR*(*fint)(x0,wl)*qe*dwl;
			}
		}
	}

    // Copy distribution to histogram
    res.hrad = new TH1D("hIR","Density of photoelectrons on radius;R, mm",Nbins,Rmin0,Rmax0);
    for(int i=0; i<Nbins; i++)
    	res.hrad->SetBinContent(i+1,I[i]);
    delete [] I;

    res.valid = kTRUE;
    res.npe = res.hrad->Integral()*dR;
    Double_t factor = (ROOT::Math::expint(res.npe)-eulergamma-log(res.npe))/(exp(res.npe)-1);
    res.R = res.hrad->GetMean();
    res.sigR = res.hrad->GetRMS();
    res.sigR_t = res.sigR*sqrt(factor);
    Double_t Lav = Ltot-La*0.5, cos2 = 1/(1+pow(res.R/Lav,2));
    res.Ang = 1000*atan(res.R/Lav); //mrad
    res.sigAng = 1000*res.sigR*cos2/Lav; //mrad
    res.sigAng_t = res.sigAng*sqrt(factor);
 
    if( draw ) {
        cout << "======Calculation results======\n" << setprecision(3)
             << "  Number of photoelectrons: " << res.npe << "\n"
             << "  Mean ring radius: " << res.R << " mm\n"
             << "  Radius single-photon sigma: " << res.sigR << " mm\n"
             << "  Radius sigma per ring: " << res.sigR_t << " mm\n"
             << "  Mean angle: " << res.Ang << " mrad\n"
             << "  Angle single-photon sigma: " << res.sigAng << " mrad\n"
             << "  Angle sigma per ring: " << res.sigAng_t << " mrad\n"
             << "===============================\n";
        res.hrad->Draw();  
    }

	return res;
}
