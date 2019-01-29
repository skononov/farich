
#include "TGraph.h"
#include "TSpline.h"
#include "TMath.h"
#include "MathFunctions.h"
#include "TGraphCorr.h"

ClassImp(TGraphCorr)

TGraphCorr::~TGraphCorr()
{}

Double_t TGraphCorr::Eval(Double_t x, TSpline* spline, Option_t* option) const
{
	if( !fNpoints ) return Nan;
	if( fNpoints==1 ) return fY[0];
	if (!spline) {
		TString opt = option;
		opt.ToLower();
		if (opt.Contains("s")) {
			// spline interpolation creating a new spline
			TSpline3 *s = new TSpline3("",this);
			Double_t result = s->Eval(x);
			delete s;
			return result;
		}
		//linear interpolation
		//find point in graph immediatly below x
		//In case x is < fX[0] or > fX[fNpoints-1] return the extrapolated point
		Int_t low = TMath::BinarySearch(fNpoints,fX,x);
		if (low == fNpoints-1) {low=fNpoints-2;}
		if (low == -1) {low=0;}
		Int_t up=low+1;
		if (fX[up]==fX[low]) return fY[low];
		Double_t yn = fY[low]*(fX[up]-x)+fY[up]*(x-fX[low]);
		return yn/(fX[up]-fX[low]);
	} else {
		//spline interpolation using the input spline
		return spline->Eval(x);
	}
}

Double_t TGraphCorr::EvalX(Double_t y, Option_t* option) const
{
	if( !fNpoints ) return Nan;
	if( fNpoints==1 ) return fY[0];
	TString opt = option;
	opt.ToLower();
	if (opt.Contains("s")) {
		// spline interpolation creating a new spline
		TSpline3 *s = new TSpline3("",fY,fX,fNpoints);
		Double_t result = s->Eval(y);
		delete s;
		return result;
	}
	//linear interpolation
	//In case y is out of range return the extrapolated point
	Int_t low;
	for(low=0; low<fNpoints-1; low++) {
		if( (y-fY[low])*(y-fY[low+1])<0 ) break;
	}
	if (low == fNpoints-1) {//intersection not found
		Double_t mindiff=TMath::Abs(y-fY[0]);
		Int_t imin=0;
		for(low=1; low<fNpoints; low++) {
			if( TMath::Abs(y-fY[low])<mindiff ) {
				mindiff=TMath::Abs(y-fY[low]);
				imin=low;
			}
		}
		return fX[imin];
	}
	Int_t up=low+1;
	if (fY[up]==fY[low]) return fX[low];
	Double_t xn = fX[low]*(fY[up]-y)+fX[up]*(y-fY[low]);
	return xn/(fY[up]-fY[low]);
}

