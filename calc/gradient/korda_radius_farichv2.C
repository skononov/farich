#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "aerdis.C"
TGraph *TQE = new TGraph("/home/skononov/work/rich/qe/hambs.dat");
TF1 *ftan = 0;

Int_t Nlambda = 50, Nbins = 200;
Double_t Lsc=50;
TH1F* hIR = new TH1F("hIR","Density of photoelectrons on radius;R, mm",Nbins,10.,30.);

Double_t func_r(Double_t *xx, Double_t *par)
{
	Double_t x0 = *xx;
	Double_t nn = par[0], k = par[1], La = par[2], Lv = par[3], b2 = par[4];
	Double_t n0 = nn+k*x0;
	ftan->SetParameters(nn,k,n0*n0);
	return ftan->Integral(x0,La) + Lv/sqrt(1/(n0*n0-1/b2)-1);
}

void derivatives(TF1* f, Double_t x, Double_t& der1, Double_t& der2)
{
	const Int_t nmaxsteps=10;
	Double_t eps = 1e-3;
	der1 = f->Derivative(x,0,eps);
	Int_t i = 0;
	while( f->DerivativeError() > 0.05*fabs(der1) ) {
		eps /= 2.;
		der1 = f->Derivative(x,0,eps);
		i++;
		if( i>=nmaxsteps ) break;
	}
	der2 = f->Derivative2(x,0,eps);
}

Double_t radius_farichv2(Double_t nn,Double_t nk, Double_t La, Double_t Ltot, Double_t b)
{
	const Double_t Const=2*TMath::Pi()/137*1e+12;
	const Double_t eps = 1e-10;

	Double_t b2 = b*b, Lv = Ltot-La;
	Double_t k = (nk-nn)/La;
	Double_t lambda_min = TQE->GetX()[0];
	Double_t lambda_max = TQE->GetX()[TQE->GetN()-1];
	Double_t dwl = (lambda_max-lambda_min)/Nlambda;

	TString stan, slam;
	stan.Form("1/sqrt(([0]+[1]*x)**2/([2]-1/%6f)-1)",b2);
	slam.Form("%6f*([0]**2-1/%6f)/[0]/x**2*exp(((400/x)**4)*[0]*%6f*([1]-%6f)/%6f)",Const,b2,b,La,Lsc);
	TF1* flam = new TF1("flam",slam,lambda_min,lambda_max);
	ftan = new TF1("ftan",stan,0,La);
	TF1* frx0 = new TF1("frx0",func_r,0.,La,5);

	frx0->SetParNames("nn","k","La","Lv","b2");
	frx0->SetParameters(nn,k,La,Lv,b2);

	Double_t Rmean = frx0->Eval(La/2);
	Int_t rbin1 = hIR->FindBin(TMath::Max(0.,Rmean-20)), rbin2 = hIR->FindBin(TMath::Min(Rmean+20,200.));
	cout << "Search radia in the range: " << hIR->GetBinCenter(rbin1) << ".." << hIR->GetBinCenter(rbin2) << endl;	 
	Double_t dR = hIR->GetBinWidth(1);

	hIR->Reset();
	hIR->Sumw2();

	for(int j=0; j<Nlambda; j++){
		Double_t wl = lambda_min+dwl*(j+0.5);
		Double_t w_nn = aerindex(nn,wl);
		Double_t w_nk = aerindex(nk,wl);
		Double_t w_k = (w_nk-w_nn)/La;

		Double_t qe = TQE->Eval(wl)/100;
		frx0->SetParameter(0,w_nn);
		frx0->SetParameter(1,w_k);
		Double_t x_ext = frx0->GetMaximumX();
		Double_t Rn = frx0->Eval(0.), Rk = frx0->Eval(La), Rmin = TMath::Min(Rn,Rk), Rmax = frx0->Eval(x_ext);
		Bool_t monotonic = kTRUE, two_roots = kFALSE;
		if( x_ext<eps || La-x_ext<eps )
			monotonic = kFALSE;

		for(int bin=rbin1; bin<=rbin2; bin++){
			Double_t R0 = hIR->GetBinCenter(bin);
			if( R0<Rmin || R0>Rmax )
				continue;

			Double_t R1 = R0+dR;
			if( R1 > Rmax )
				R1 = Rmax;

			Double_t x1=-1, x2=-1;
			if( monotonic ) {			
				x1 = frx0->GetX(R0);
				if( x1<eps || La-x1<eps )
					continue;
			} else {
				x1 = frx0->GetX(R0,0.,x_ext);
				x2 = frx0->GetX(R0,x_ext,La);
				if( x1<eps || La-x1<eps ) {
					if( x2<eps || La-x2<eps )
						continue;
					x1 = x2;
				} else if( x2>=eps && La-x2>=eps ) {
					two_roots = kTRUE;
				}
			}

			Double_t I = 0., dRdx, d2Rdx;
			derivatives(frx0,x1,dRdx,d2Rdx);
			Double_t D = dRdx*dRdx+2*d2Rdx*(R1-R0);
			if( D < 0 ) 
				continue;
			Double_t dxdR = fabs((sqrt(D)-fabs(dRdx))/d2Rdx/dR);
			Double_t n = w_nn + w_k*x1;
			flam->SetParameters(n,x1);
			I = dxdR*(*flam)(wl)*qe*dwl;
			if( two_roots ) {
				derivatives(frx0,x2,dRdx,d2Rdx);
				D = dRdx*dRdx+2*d2Rdx*(R1-R0);
				if( D < 0 ) 
					continue;
				dxdR = fabs((sqrt(D)-fabs(dRdx))/d2Rdx/dR);
				n = w_nn + w_k*x2;
				flam->SetParameters(n,x2);
				I += dxdR*(*flam)(wl)*qe*dwl;
			}

			hIR->SetBinContent(bin,hIR->GetBinContent(bin)+I);
		}
	}

	hIR->Draw();  

	return hIR->GetRMS();
}
