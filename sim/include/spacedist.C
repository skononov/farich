#include "TF1.h"

/**************************************************************
The function for resolution studies - convolution of
linear step and gaussian functions plus constant background
  Has 6 free parameters:
  1: Signal level in histogram
  2: Width of the distribution
  3: Linear coefficient
  4: position of the distribution's center
  5: Coordinate resolution - sigma of gaussian
  6: Background level in histogram
***************************************************************/

Double_t spacedist(Double_t *xx,Double_t *par)
{
	static const Double_t sqrt2PI=sqrt(2*TMath::Pi()), sqrt2=TMath::Sqrt(2.);
	Double_t C=par[0], w=par[1], k=par[2], m=par[3], sig=par[4], B=par[5];
	Double_t x=*xx, arg1=(x-w/2-m)/sqrt2/sig, arg2=(x+w/2-m)/sqrt2/sig;
	return B + C*( 0.5*(k*(x-m)+1/w)*(TMath::Erf(arg2)-TMath::Erf(arg1)) -
		           k*sig/sqrt2PI*(exp(-arg1*arg1)-exp(-arg2*arg2)) );
}

TF1* fsp1=new TF1("fsp1",spacedist,0,1,6);

/**************************************************************
 Has 8 free parameters:
  1: Constant
  2: Width of the distribution
  3: Cubic coefficient
  4: Quadratic coefficient
  5: Linear coefficient
  6: Position of the distribution's center
  7: Sigma of gaussian 1
  8: Sigma of gaussian 2
***************************************************************/

Double_t spacedist3(Double_t *xx,Double_t *par)
{
	static const Double_t sqrt2=TMath::Sqrt(2.);
	Double_t x=*xx;
	Double_t C=par[0], hw=par[1]/2, a=par[2], b=par[3], c=par[4],
		     m=par[5], sig1=par[6], sig2=par[7];
	Double_t t=x-m;
	Double_t pol3=a*t*t*t+b*t*t+c*t+1;

	if( C==0.0 ) return 0;

	return C*pol3*pol3*(1+TMath::Erf((t+hw)/sqrt2/sig1))*
        (1-TMath::Erf((t-hw)/sqrt2/sig2));
}

TF1* fsp3=new TF1("fsp3",spacedist3,0,1,8);

/**************************************************************
 Has 8 free parameters:
  1: Constant
  2: Width of the distribution
  3: Cubic coefficient
  4: Quadratic coefficient
  5: Linear coefficient
  6: Free coefficient
  7: Position of the distribution's center
  8: Sigma of gaussian 1
  9: Sigma of gaussian 2
***************************************************************/

Double_t spacedist3at(Double_t *xx,Double_t *par)
{
	static const Double_t sqrt2=TMath::Sqrt(2.), halfPi=TMath::PiOver2();
	Double_t x=*xx;
	Double_t C=par[0], hw=par[1]/2, a=par[2], b=par[3], c=par[4], d=par[5],
		     m=par[6], sig1=par[7], sig2=par[8];
	Double_t t=x-m;

	if( C==0.0 ) return 0;

/*  old formula
	return C*(halfPi+atan(a*t*t*t+b*t*t+c*t+d))/(halfPi+atan(d))*
		(1+TMath::Erf((t+hw)/sqrt2/sig1))*
		(1-TMath::Erf((t-hw)/sqrt2/sig2));*/

	Double_t atNorm=2*halfPi*hw;
	Double_t arg1=c*hw+d, arg2=-c*hw+d;
	if( c!=0.0 )
		atNorm+=(arg1*atan(arg1)-arg2*atan(arg2)-0.5*log((1+arg1*arg1)/(1+arg2*arg2)))/c;
	else
		atNorm+=2*atan(d)*hw;

	return C*(halfPi+atan(c*t+d))/atNorm*
		(1+TMath::Erf((t+hw)/sqrt2/sig1))*
        (1-TMath::Erf((t-hw)/sqrt2/sig2));
}

TF1* fsp3at=new TF1("fsp3at",spacedist3at,0,1,9);


