
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TKey.h"
#include "TMath.h"
#include "TSpline.h"
#include "TGraphCorr.h"

template class vector<TGraphCorr*>;
template class vector<TSpline3*>;

Int_t Nr=0, Nphi=0;
Double_t dPhi=2*TMath::Pi();

vector<TGraphCorr*> vgBeta, vgSigBeta, vgRingBound;
vector<TSpline3*> vsBeta;

TDatime modDatime(2010,4,30,13,0,0);
Bool_t newCal=kTRUE;

void graphs_clear()
{
	size_t i;
	for(i=0; i<vgBeta.size(); i++)
		delete vgBeta[i];
	for(i=0; i<vgSigBeta.size(); i++)
		delete vgSigBeta[i];
	for(i=0; i<vgRingBound.size(); i++)
		delete vgRingBound[i];
	vgBeta.clear();
	vgSigBeta.clear();
	vgRingBound.clear();
	vsBeta.clear();
}

TGraph *loadgraph(const char* name)
{
	TGraph* g = (TGraph*)gDirectory->FindObject(name);
	if (!g) {
		g = (TGraph*)gDirectory->Get(name);
		gDirectory->Append(g);
	}
	return g;
}

Int_t loadcal(const char* calfn)
{
	TFile *calfile = new TFile(calfn);
	if (!calfile->IsOpen()) return 0;

	std::cout<<"Reading beta calibrations from "<<calfile->GetName()<<std::endl;
	calfile->cd();

	TGraph *aGraph=0;
	TGraphCorr *cGraph=0;

	const char *bfmt, *sbfmt, *boundfmt;

	char name[20];

	Nphi=0;
	Nr=0;
	TDatime cdt;

	if (loadgraph("beta1_0")) { //angle mode
		cdt=gDirectory->GetKey("beta1_0")->GetDatime();
		TGraph *gr, *gphi;
		while(1) { //get number of rings and phi divisions
			sprintf(name,"beta%d_0",Nr+1);
			if ((gr=loadgraph(name))) Nr++;
			sprintf(name,"beta1_%d",Nphi);
			if ((gphi=loadgraph(name))) Nphi++;
			if (!gr && !gphi) break;
		}
		bfmt="beta%d_%d";
		sbfmt="sigbeta%d_%d";
		boundfmt="bound%d_%d";
	} else {
		Nphi=1;
		while(1) { //get number of rings
			sprintf(name,"beta%d",Nr+1);
			if (!loadgraph(name)) break;
			Nr++;
		}
		bfmt="beta%d";
		sbfmt="sigbeta%d";
		boundfmt="bound%d";
	}
	if (Nr==0) {
		cerr<<"Calibration graphs not found"<<std::endl;
		delete calfile;
		return 0;
	}

	if (Nr==1)  std::cout<<"Single ring mode determined.";
	else        std::cout<<"Multi-ring mode: "<<Nr<<" rings.";
	if (Nphi>1) std::cout<<" "<<Nphi<<" phi divisions.";
	std::cout<<std::endl;

	graphs_clear();
	vgBeta.reserve(Nr*Nphi);
	vgSigBeta.reserve(Nr*Nphi);
	vgRingBound.reserve((Nr-1)*Nphi);
	vsBeta.reserve(Nr*Nphi);
	for (Int_t ir=0; ir<Nr; ir++) {
		for (Int_t iphi=0; iphi<Nphi; iphi++) {
			sprintf(name,bfmt,ir+1,iphi);
			aGraph = loadgraph(name);
			if (!aGraph) break;
			cGraph = new TGraphCorr(*aGraph);
			vgBeta.push_back(cGraph);
			vsBeta.push_back(new TSpline3("",aGraph));

			if (Nr==1 && Nphi==1) break;

			sprintf(name,sbfmt,ir+1,iphi);
			aGraph = loadgraph(name);
			if (!aGraph) break;
			cGraph = new TGraphCorr(*aGraph);
			vgSigBeta.push_back(cGraph);

			if (Nr==1 || ir>Nr-2) continue;

			sprintf(name,boundfmt,ir+1,iphi);
			aGraph = loadgraph(name);
			if (!aGraph) break;
			cGraph = new TGraphCorr(*aGraph);
			vgRingBound.push_back(cGraph);
		}
		if (!aGraph) break;
	}

	delete calfile;

	if (!aGraph) {
		cerr<<"Graph with name "<<name<<" not found in calibration file"<<std::endl;
		return 0;
	}

	if (Nphi>1) {
		dPhi=TMath::Pi()/Nphi;
		if( cdt.Convert() > modDatime.Convert() ) {
			cout<<"New phi-stepping for calibration file after "<<modDatime.AsString()<<endl;
			dPhi*=2;
		} else {
			newCal=kFALSE;
			cout<<"Old phi-stepping for calibration file prior to "<<modDatime.AsString()<<endl;
		}
	}

	return Nr*Nphi;
}

static Double_t radius_to_beta(Double_t *x,Double_t *par)
{
	Int_t ir=int(par[0]); //ring index
	Int_t iphi=int(par[1]); //phi division index

	if (ir<0 || ir>=Nr || iphi<0 || iphi>=Nphi) return 0.0;

	Double_t &r=*x; //radius value
	return vgBeta[ir*Nphi+iphi]->Eval(r);
}

static Double_t sigma_beta(Double_t *x,Double_t *par)
{
	Int_t ir=int(par[0]); //ring index
	Int_t iphi=int(par[1]); //phi division index

	if (ir<0 || ir>=Nr || iphi<0 || iphi>=Nphi) return 0.0;

	Double_t &beta=*x;
	return vgSigBeta[ir*Nphi+iphi]->Eval(beta);
}

static Double_t ring_number(Double_t *x,Double_t *par)
{
	if (Nr==1) return 0;

	Double_t beta=par[0]; //velocity of the particle
	Double_t &r=*x;

	Int_t ir;
	for (ir=0; ir<Nr-1; ir++) {
		if (r<vgRingBound[ir]->EvalX(beta)) break;
	}

	return ir;
}

