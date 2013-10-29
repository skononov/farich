#include <cmath>
#include "Minuit2FCN.h"
#include "MLADescription.h"
#include "Spectrum.h"

extern MLADescription mla;
extern Spectrum phdeteff;
extern float pixelsize;

double Minuit2FCN::operator()(const std::vector<double>& x) const
{
	double Tnorm=0;
	for(int l=0; l<mla.GetNlayers(); l++)
		Tnorm+=x[2*l];
	double T=mla.GetTotalThickness();
	for(int l=0; l<mla.GetNlayers(); l++) {
		if( l>0 ) mla.SetIndex(l,x[2*l-1]);
		mla.SetThickness(l,x[2*l]*T/Tnorm);
	}

	MLAResult& res=mla.Calculate(phdeteff,pixelsize,mla.GetBeta());

	if( !res.valid ) return NAN;

	return res.sigma_t;
}

