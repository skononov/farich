
#include <vector>
#include <iostream>
#include "TMultiGraph.h"
#include "TFile.h"
#include "TMath.h"

#include "TGraphCorr.h"
#include "spacedist.C"

using namespace std;

Int_t Nr=0, Np=0;
Double_t ***vvPars=0;
Double_t* vBeta=0;

static void freePars()
{
	if( vvPars ) {
		for(Int_t l=0; l<Nr; l++) {
			for(Int_t p=0; p<Np; p++)
				delete [] vvPars[l][p];
			delete [] vvPars[l];
		}
        delete [] vvPars;
	}
	vvPars=0;
	if( vBeta ) delete vBeta;
	vBeta=0;
}

static void allocatePars()
{
	vvPars=new Double_t**[Nr];
	for(Int_t l=0; l<Nr; l++) {
		vvPars[l]=new Double_t*[Np];
		for(Int_t p=0; p<Np; p++)
			vvPars[l][p]=new Double_t[9];
	}
	vBeta=new Double_t[Np];
}


Int_t loadcal(const char* calfn)
{
	freePars();

	TFile calfile(calfn);
	if (!calfile.IsOpen()) return 0;

	cout<<"Reading pdf calibrations from "<<calfile.GetName()<<endl;
	calfile.cd();

	TMultiGraph* mgRaw[9];

	char name[20];

	for(Int_t ipar=0; ipar<9; ipar++) {
		sprintf(name,"mgpar%d",ipar);
		mgRaw[ipar]=(TMultiGraph*)calfile.Get(name);
		if( !mgRaw[ipar] ) {
			cerr<<"TMultiGraph object for parameter "<<ipar<<" missing"<<endl;
			return 0;
		}
		if( ipar==0 ) {
			Nr=mgRaw[0]->GetListOfGraphs()->GetSize();
		} else if( Nr!=mgRaw[ipar]->GetListOfGraphs()->GetSize() ) {
			cerr<<"TMultiGraph for parameter "<<ipar<<" is not compatible with number of rings defined"<<endl;
			return 0;
		}
	}

	//Determine points on beta to add, include beta=1 point
	Np=0;
	vector<Double_t> vnewbeta;
	Double_t **vmean=new Double_t*[Nr];
	Double_t **vwidth=new Double_t*[Nr];
	for (Int_t l=0; l<Nr; l++) {
		vmean[l]=((TGraphCorr*)mgRaw[6]->GetListOfGraphs()->At(l))->GetY();
		vwidth[l]=((TGraphCorr*)mgRaw[1]->GetListOfGraphs()->At(l))->GetY();
	}

	TGraphCorr* g=(TGraphCorr*)mgRaw[0]->GetListOfGraphs()->At(0);
	Int_t Np_old=g->GetN();
	Double_t *vbeta=g->GetX();

	for(Int_t p=0; p<Np_old-1; p++) {
		vnewbeta.push_back(vbeta[p]);
		Int_t ninrange=0;
		for (Int_t l=0; l<Nr; l++) {
			if( vmean[l][p]==0 ) continue;
			Int_t n=(Int_t)(1+10*(vmean[l][p+1]-vmean[l][p])/(vwidth[l][p+1]+vwidth[l][p]));
			if( n>ninrange ) ninrange=n;
		}
		if( !ninrange ) continue;
		for(Int_t i=1; i<ninrange; i++) {
            vnewbeta.push_back((vbeta[p+1]*i+vbeta[p]*(ninrange-i))/ninrange);
		}
	}
	vnewbeta.push_back(vbeta[Np_old-1]);
	vnewbeta.push_back(1.0);

	delete [] vmean;
	delete [] vwidth;

	Np=vnewbeta.size();

	cout<<Nr<<" rings "<<Np<<" points in calibration"<<endl;

	allocatePars();

	for(Int_t p=0; p<Np; p++) vBeta[p]=vnewbeta[p];

	//Make finer calibration graphs in beta
	for (Int_t ipar=0; ipar<9; ipar++) {
		for (Int_t l=0; l<Nr; l++) {
			TGraphCorr* gpar=(TGraphCorr*)mgRaw[ipar]->GetListOfGraphs()->At(l);
			Double_t* vpar0=((TGraphCorr*)mgRaw[0]->GetListOfGraphs()->At(l))->GetY();

			Int_t p_old=0;
			for(Int_t p=0; p<Np; p++) {
				if( vbeta[p_old+1]<=vnewbeta[p] ) p_old++;
				if( vpar0[p_old]==0 )
					vvPars[l][p][ipar]=0;
				else
					vvPars[l][p][ipar]=gpar->Eval(vnewbeta[p]);
			}
		}
	}

	for (Int_t ipar=0; ipar<9; ipar++) delete mgRaw[ipar];

	calfile.Close();

	//Normalize interpolated functions
	for(Int_t p=0; p<Np; p++) {
		Double_t integral=0;
		for(Int_t l=0; l<Nr; l++) {
			Double_t *par=vvPars[l][p];
			fsp3at->SetParameters(par);
			Double_t x1=par[6]-par[1]/2-5*par[7], x2=par[6]+par[1]/2+5*par[8];
			integral+=fsp3at->Integral(x1,x2);
		}
		if( integral>0 ) {
			for(Int_t l=0; l<Nr; l++)
				vvPars[l][p][0]/=integral;
		}
	}

	return 1;
}

Double_t pdfvalp(Int_t p,Double_t r,Int_t layer=-1)
{
	Double_t val=0;
	if( layer>=0 )
		val=spacedist3at(&r,vvPars[layer][p]);
	else {
		for(Int_t l=0; l<Nr; l++)
			val+=spacedist3at(&r,vvPars[l][p]);
	}

	return val;
}

static Double_t betaLast=0;
static Int_t pointLast=-1;

Double_t pdfval(Double_t *x,Double_t *par)
{
	Double_t &beta=par[0];
	Int_t layer=(Int_t)par[1];

	Int_t p;
	if( betaLast==beta )
		p=pointLast;
	else {
		p=TMath::BinarySearch(Np,vBeta,beta);
		if( p<0 )
			p=0;
		else if( p<Np-1 )
			if( vBeta[p+1]-beta<beta-vBeta[p] )
				p++;
		pointLast=p;
	}

	return pdfvalp(p,*x,layer);
}


Double_t findMaxLikelihood(Int_t n,Double_t *r,Int_t layer=-1)
{
	Double_t lhmax=0, lhmax_p=0, lhmax_n=0, lh=0, lhprev=0;
	Int_t popt=-1;

	for(Int_t p=0; p<Np; p++) {
		lhprev=lh;
		lh=1;

		for(Int_t i=0; i<n; i++)
			lh*=pdfvalp(p,r[i],layer);

		if( lhmax<lh ) {
			lhmax_p=lhprev;
			lhmax=lh;
			popt=p;
		} else {
			if( popt==p-1 ) lhmax_n=lh;
		}
	}

	Double_t beta=vBeta[popt];

	if( popt==0 || popt==Np-1 ) return beta;

	//make quadratic approximation
	Double_t t1=vBeta[popt-1]-vBeta[popt], t2=vBeta[popt+1]-vBeta[popt];
	Double_t y1=lhmax_p-lhmax, y2=lhmax_n-lhmax;

	beta+=0.5*(y1*t2*t2-y2*t1*t1)/(y1*t2-y2*t1);

	return beta;
}

Double_t findOptimalBeta(Double_t r,Double_t beta_est,Int_t &layer)
{
	Double_t pdfmax=0, dmin=1e3;
	Int_t lpdf=-1, ldist=-1;
	layer=-1;

	Double_t par[2]={beta_est,-1};

	for(Int_t l=0; l<Nr; l++) {
		par[1]=l;
		Double_t pdf=pdfval(&r,par);
		Int_t p=pointLast;
		if( vvPars[l][p][0]>0 ) {
			Double_t dist=r-vvPars[l][p][6];
			Double_t d=fabs(dist/(vvPars[l][p][1]/2+(dist<0?vvPars[l][p][7]:vvPars[l][p][8])));
			if( d<dmin ) {
				dmin=d;
				ldist=l;
			}
		}
		if( pdfmax<pdf ) {
			pdfmax=pdf;
			lpdf=l;
		}
	}
	if( lpdf>=0 ) layer=lpdf;
	else if( ldist>=0 ) layer=ldist;

	return findMaxLikelihood(1,&r,layer);
}


