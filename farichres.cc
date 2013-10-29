#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <libgen.h>
#include <string>

#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TMath.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include "MLADescription.h"
#include "Spectrum.h"
#include "Minuit2FCN.h"

using namespace std;

static const char *optstring="qn:s:N:D:T:p:b:e:o:O:m";
static const char *progname;

//��������� ���������
static const char *qefn;
static bool  batch=false;
static float ri1=1.07;
static float Lsc=50.;
static int   nlayers=3;
static float D=100.;
static float T=25.;
static float efficiency=1.0;
float pixelsize=0;
static float opbeta=1.0;
static bool  minimize=false;
static const char *outfn="";
static const char *macfn="";

MLADescription mla;
Spectrum phdeteff;

double get_max_sensitivity_wl()
{
	double wl1, wl2;
	phdeteff.GetRange(wl1,wl2);
	double wlstep=(wl2-wl1)/50;
	double wl0=400, smax=0;

	for(int i=0; i<=50; i++) {
		double wl=wl1+i*wlstep;
		double rn=MLADescription::AerogelRefIndex(ri1,400.,wl)*opbeta;
		if( rn<=1.0 ) continue;
		double Latt=Lsc*pow(wl/400,4);
		double s=(1-1/(rn*rn))*Latt*phdeteff.Evaluate(wl)*(1-exp(-T*rn/Latt))/rn;
		if( smax<s ) {
			smax=s;
			wl0=wl;
		}
	}

	return wl0;
}

void Usage(int status)
{
	cout<<"Usage: "<<progname<<" [OPTIONS] qefile\n"
		<<"��������� ���������� ��������� ������������ ����� � ������������ ������������ "
		<<"���������. ������������ ����������. �������� qefile ������ ���� � ����� � "
		<<"������� � ��������� ������������� ��������� ���������\n"
		<<" OPTIONS:\n"
		<<"   -q            �������� ����� ��� ������� � ��������������\n"
		<<"   -e eff        ������ � ������������� ��������� ��������� 0<eff<=1 ("<<efficiency<<")\n"
		<<"   -p size       ������ ����������� ������� ��������� ���������, �� ("<<pixelsize<<")\n"
		<<"   -D distance   ���������� �� ������ ��������� �� ��������� ���������, �� ("<<D<<")\n"
		<<"   -n ri1        ������������ ���������� ����������� ��������� �� 400 �� ("<<ri1<<")\n"
		<<"   -N nlayers    ����� ����� �������� ("<<nlayers<<")\n"
		<<"   -T thickness  ������� ���������, �� ("<<T<<")\n"
		<<"   -s Lsc        ����� ��������� �� 400 ��, �� ("<<Lsc<<")\n"
		<<"   -b beta       �������������� �������� ��� ������ �������� ("<<opbeta<<")\n"
		<<"   -o filename   �������� root ���� ��� ���������� ������������� �� ������� ("<<outfn<<")\n"
		<<"   -O filename   ��������� �������� ��������� � ��������� ���� Geant4 ("<<macfn<<")\n"
		<<"   -m            ����������� ������ �� ������� �� ����\n"
		<<endl;
	exit(status);
}

int main(int argc, char* argv[])
{
	progname=argv[0];

	if( argc==1 ) Usage(0);

//=========��������� ���������� ���������==========//
	int opt;
	while( (opt=getopt(argc,argv,optstring)) > 0 ) {
		if( opt=='?' )
			Usage(1);
		else if( opt=='q' ) {
			batch=true;
		}
		else if( opt=='n' ) {
			float n=atof(optarg);
			if( n<1 || n>1.41 ) {
				cerr<<optarg<<": ������������ ���������� �����������"<<endl;
				return 1;
			}
			ri1=n;
		}
		else if( opt=='s' ) {
			float l=atof(optarg);
			if( l<=0 ) {
				cerr<<optarg<<": ����������� ������ ����� ��������"<<endl;
				return 1;
			}
			Lsc=l;
		}
		else if( opt=='N' ) {
			int n=atoi(optarg);
			if( n<1 ) {
				cerr<<optarg<<": ����������� ������ ����� �����"<<endl;
				return 1;
			}
			nlayers=n;
		}
		else if( opt=='D' ) {
			float d=atof(optarg);
			if( d<=0.0 ) {
				cerr<<optarg<<": ����������� ������ ���������� �� ��������� �� ���������"<<endl;
				return 1;
			}
			D=d;
		}
		else if( opt=='T' ) {
			float t=atof(optarg);
			if( t<=0.0 ) {
				cerr<<optarg<<": ����������� ������ ������� ���������"<<endl;
				return 1;
			}
			T=t;
		}
		else if( opt=='b' ) {
			float b=atof(optarg);
			if( b<=0.0 || b>1.0 ) {
				cerr<<optarg<<": ����������� ������ �������� ������� ��� ����������� ���������"<<endl;
				return 1;
			}
			opbeta=b;
		}
		else if( opt=='e' ) {
			float e=atof(optarg);
			if( e<=0.0 || e>1.0 ) {
				cerr<<optarg<<": ����������� ����� ������ ������������� ��������� ���������"<<endl;
				return 1;
			}
			efficiency=e;
		}
		else if( opt=='p' ) {
			float s=atof(optarg);
			if( s<=0.0 ) {
				cerr<<optarg<<": ����������� ����� ������ �������"<<endl;
				return 1;
			}
			pixelsize=s;
		}
		else if( opt=='m' ) {
			minimize=true;
		}
		else if( opt=='o' ) {
			outfn=optarg;
		}
		else if( opt=='O' ) {
			macfn=optarg;
		}
		else {
			cerr<<opt<<": ����������� �����"<<endl;
			Usage(1);
		}
	}
	if( T >= D ) {
		cerr<<"T="<<T<<", D="<<D
			<<": ������� ��������� ������ ���������� �� ��������� ���������"<<endl;
		return 1;
	}
	if( ri1*opbeta < 1.0 ) {
		cerr<<"ri1="<<ri1<<" beta="<<opbeta
			<<": ����������� �������� � ������ ����"<<endl;
		return 1;
	}
	if( optind >= argc )
		Usage(1);

	qefn=argv[optind];

//==========������========//
	cout<<"________________________________________________\n"
		<<"���������� ���������� ������������� ���������.\n"
		<<"���������:\n"
		<<"  ���� � ��������� ��������������:\n   "<<qefn<<"\n"
		<<"  ������ ������������� ���:                     "<<efficiency<<"\n"
		<<"  ������ ������� ���:                           "<<pixelsize<<" ��\n"
		<<"  ���������� �� ��������� �� ���:               "<<D<<" ��\n"
		<<"  ���������� ����������� ������� ���� �� 400��: "<<ri1<<"\n"
		<<"  ����� ����� ���������:                        "<<nlayers<<"\n"
		<<"  ������ ������� ���������:                     "<<T<<" ��\n"
		<<"  ����� ��������� � �������� �� 400 ��:         "<<Lsc<<" ��\n"
		<<"  ����������� ��� ��������:                     "<<opbeta<<"\n"
		<<"  ����������� "<<(minimize?"��������":"���������")<<"\n";
	if( *outfn ) cout<<"  �������� root-����: "<<outfn<<"\n";
	if( *macfn ) cout<<"  ���� ������� ��� Geant4: "<<macfn<<"\n";
	cout<<"________________________________________________"<<endl;

	TRint *app=0;
	if( !batch ) {
		int ac=2;
		char *av[2]={"rint","-l"};
		app=new TRint("rint",&ac,av);
	}

	phdeteff.ReadFile(qefn);
	if( phdeteff.IsEmpty() )
		return 1;

	if( efficiency<1.0 ) phdeteff.Scale(efficiency);

	double wl1, wl2;
	phdeteff.GetRange(wl1,wl2);
	double wl0=get_max_sensitivity_wl();

	cout.precision(3);
	cout<<"��������� ������������� ���������� �� "<<wl1<<" �� "<<wl2<<" ��, ����� "
		<<phdeteff.Size()<<" ��������."<<endl;
	cout<<"����� ����� ������������ ���������������� � ��: "<<wl0<<" ��"<<endl;
	cout.precision(0);

	//������� ������������ ����������� �������� �� ��������� ������ ��� ��������� ���
	//����� ����� ������������ ����������������
	double t0=D-T; //proximity distance
	MLADescription mla0(t0,opbeta,wl0);
	mla0.SetScatteringLength(Lsc);

	mla0.MakeFixed(nlayers,D,ri1);
//	for(int l=0; l<nlayers; l++)
//		mla0.AddAlayer(ri1,T/nlayers);

	cout<<"�������� ����������� ��������:"<<endl;
	mla0.Print("  ");

	MLAResult& res=mla0.Calculate(phdeteff,pixelsize);

	cout<<"  ������� ����� ��������������: "<<res.npe<<endl;
	cout<<"  ������� ������:               "<<res.radius<<endl;
	cout<<"  ������ ������� �� 1 �����:    "<<res.sigma1<<endl;
	cout<<"  ������ ������� �� ����:       "<<res.sigma_t<<endl;

	mla=mla0;
	if( minimize ) {
		Minuit2FCN fcn(0.1*res.sigma_t);

		char name[20];
		ROOT::Minuit2::MnUserParameters pars;
		pars.Add(string("t1"),mla.GetThickness(0),0.01,0.0,T);
		for(int l=1; l<nlayers; l++) {
			sprintf(name,"n%d",l+1);
			pars.Add(string(name),mla.GetIndex(l),0.005,1.0,2.0);
			sprintf(name,"t%d",l+1);
			pars.Add(string(name),mla.GetThickness(l),0.5,0.0,T);
		}

		ROOT::Minuit2::MnStrategy strategy(1);

		ROOT::Minuit2::MnMinimize combmin(fcn, pars, strategy);

//		combmin.SetPrecision(1e-8);

		cout<<"������ �������"<<endl;

		ROOT::Minuit2::FunctionMinimum minimum = combmin(1000,0.1);

		cout<<minimum;

		ROOT::Minuit2::MnUserParameterState state=minimum.UserState();

		if( !minimum.IsValid() ) {
			cout<<"�� ������� �������������� 1 :("<<endl;
			return 0;
		}

		cout<<"����������� ���������"<<endl;
		state.SetError((unsigned int)0,0.05);
		state.RemoveLimits((unsigned int)0);
		for(int l=1; l<nlayers; l++) {
			state.SetError((unsigned int)2*l,0.05);
			state.RemoveLimits((unsigned int)2*l);
			state.SetError((unsigned int)2*l-1,0.001);
			state.RemoveLimits((unsigned int)2*l-1);
		}

		strategy.SetLowStrategy();
		ROOT::Minuit2::MnMigrad migrad(fcn, state, strategy);

		minimum=migrad(1000,0.1);

		cout<<minimum;

		if( !minimum.IsValid() ) {
			cout<<"�� ������� �������������� 2 :("<<endl;
		} else {
			state=minimum.UserState();
		}

		double Tnorm=0;
		for(int l=0; l<mla.GetNlayers(); l++)
			Tnorm+=state.Value(2*l);

		for(int l=0; l<mla.GetNlayers(); l++) {
			if( l>0 ) mla.SetIndex(l,state.Value(2*l-1));
			mla.SetThickness(l,state.Value(2*l)*T/Tnorm);
		}

		cout<<"����������� �������� ����������������"<<endl;
		mla.Print("  ");

		res=mla.Calculate(phdeteff,pixelsize);
		cout<<"  ������� ����� ��������������: "<<res.npe<<endl;
		cout<<"  ������� ������:               "<<res.radius<<endl;
		cout<<"  ������ ������� �� 1 �����:    "<<res.sigma1<<endl;
		cout<<"  ������ ������� �� ����:       "<<res.sigma_t<<endl;
	}

	TH1D hrad("hrad","Radius photoelectron distribution;radius, mm",Nr,res.rmin,res.rmax);
	for(int i=0; i<Nr; i++) hrad.SetBinContent(i+1,res.s[i]);

	if( *outfn ) {
		cout<<"��������� ������������� ������� �� ������� � "<<outfn<<endl;
		TFile *outfile=new TFile(outfn,"RECREATE");
		if( outfile->IsOpen() )
			hrad.Write("hrad");
		else
			cerr<<outfn<<": �� ���� ������� ���� �� ������"<<endl;
		outfile->Close();
		delete outfile;
	}

	if( *macfn ) {
		cout<<"��������� �������� RICH � ��������� ���� Geant4 "<<macfn<<endl;
		ofstream macfile(macfn);
		if( macfile.is_open() ) {
			macfile<<"/Rich/pmt/qeDataFile "<<qefn<<"\n"
				<<"/Rich/pmt/detection "<<efficiency<<"\n"
				<<"/Rich/radiator/scatterLength "<<Lsc<<" mm\n"
				<<"/Rich/radiator/absDataFile none\n"
				<<"/Rich/mode manual\n"
				<<"/Rich/proximity "<<t0<<" mm\n";
			for(int l=0; l<nlayers; l++) {
				macfile<<"/Rich/radiator/addLayer "<<setprecision(5)<<mla.GetIndex(l)
					<<" "<<mla.GetThickness(l)<<" mm\n";
			}
			macfile<<"/Rich/update"<<endl;
		} else {
			cerr<<macfn<<": �� ���� ������� ���� �� ������"<<endl;
		}
	}

	if( !batch ) {
		TCanvas *c=new TCanvas("c1","������������� �������������� �� �������");

		hrad.Draw("l");

		app->Run(kTRUE);
	}

	return 0;
}


