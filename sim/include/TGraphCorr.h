
#ifndef TGraphCorr_h
#define TGraphCorr_h

#include "TGraph.h"

class TGraphCorr : public TGraph {
public:
	TGraphCorr() : TGraph() {}
	TGraphCorr(Int_t n) : TGraph(n) {}
	TGraphCorr(Int_t n, const Int_t* x, const Int_t* y) :
		TGraph(n,x,y)
	{}
	TGraphCorr(Int_t n, const Float_t* x, const Float_t* y) :
		TGraph(n,x,y)
	{}
	TGraphCorr(Int_t n, const Double_t* x, const Double_t* y) :
		TGraph(n,x,y)
	{}
	TGraphCorr(const char* fn, const char* fmt = "%lg %lg", Option_t* opt = "") :
		TGraph(fn,fmt,opt)
	{}
	TGraphCorr(const TGraph& gr) :
		TGraph(gr)
	{}
	TGraphCorr(const TGraphCorr& gr) :
		TGraph(gr)
	{}
	~TGraphCorr();

	Double_t Eval(Double_t x, TSpline* spline = 0, Option_t* option="") const;
	Double_t EvalX(Double_t y, Option_t* option="") const;

	ClassDef(TGraphCorr,2);
};
#endif // TGraphCorr_h

