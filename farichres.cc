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

//Параметры программы
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
		<<"Вычисляет разрешение детектора черенковских колец с многослойным фокусирующим "
		<<"аэрогелем. Оптимизирует разрешение. Параметр qefile задает путь к файлу с "
		<<"данными о квантовой эффективности фотонного детектора\n"
		<<" OPTIONS:\n"
		<<"   -q            Задачный режим без графики и интерпретатора\n"
		<<"   -e eff        Фактор к эффективности фотонного детектора 0<eff<=1 ("<<efficiency<<")\n"
		<<"   -p size       Размер квадратного пикселя фотонного детектора, мм ("<<pixelsize<<")\n"
		<<"   -D distance   Расстояние от начала радиатора до фотонного детектора, мм ("<<D<<")\n"
		<<"   -n ri1        Максимальный показатель преломления радиатора на 400 нм ("<<ri1<<")\n"
		<<"   -N nlayers    Число слоев аэрогеля ("<<nlayers<<")\n"
		<<"   -T thickness  Толщина радиатора, мм ("<<T<<")\n"
		<<"   -s Lsc        Длина рассеяния на 400 нм, мм ("<<Lsc<<")\n"
		<<"   -b beta       Оптимизировать радиатор для данной скорости ("<<opbeta<<")\n"
		<<"   -o filename   Выходной root файл для сохранения распределения по радиусу ("<<outfn<<")\n"
		<<"   -O filename   Сохранить описание детектора в командный файл Geant4 ("<<macfn<<")\n"
		<<"   -m            Минимизация ошибки по радиусу на трек\n"
		<<endl;
	exit(status);
}

int main(int argc, char* argv[])
{
	progname=argv[0];

	if( argc==1 ) Usage(0);

//=========Обработка параметров программы==========//
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
				cerr<<optarg<<": сомнительный показатель преломления"<<endl;
				return 1;
			}
			ri1=n;
		}
		else if( opt=='s' ) {
			float l=atof(optarg);
			if( l<=0 ) {
				cerr<<optarg<<": неправильно задана длина рассения"<<endl;
				return 1;
			}
			Lsc=l;
		}
		else if( opt=='N' ) {
			int n=atoi(optarg);
			if( n<1 ) {
				cerr<<optarg<<": неправильно задано число слоев"<<endl;
				return 1;
			}
			nlayers=n;
		}
		else if( opt=='D' ) {
			float d=atof(optarg);
			if( d<=0.0 ) {
				cerr<<optarg<<": неправильно задано расстояние от радиатора до детектора"<<endl;
				return 1;
			}
			D=d;
		}
		else if( opt=='T' ) {
			float t=atof(optarg);
			if( t<=0.0 ) {
				cerr<<optarg<<": неправильно задана толщина радиатора"<<endl;
				return 1;
			}
			T=t;
		}
		else if( opt=='b' ) {
			float b=atof(optarg);
			if( b<=0.0 || b>1.0 ) {
				cerr<<optarg<<": неправильно задана скорость частицы для оптимизации радиатора"<<endl;
				return 1;
			}
			opbeta=b;
		}
		else if( opt=='e' ) {
			float e=atof(optarg);
			if( e<=0.0 || e>1.0 ) {
				cerr<<optarg<<": неправильно задан фактор эффективности фотонного детектора"<<endl;
				return 1;
			}
			efficiency=e;
		}
		else if( opt=='p' ) {
			float s=atof(optarg);
			if( s<=0.0 ) {
				cerr<<optarg<<": неправильно задан размер пикселя"<<endl;
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
			cerr<<opt<<": неизвестная опция"<<endl;
			Usage(1);
		}
	}
	if( T >= D ) {
		cerr<<"T="<<T<<", D="<<D
			<<": толщина радиатора больше расстояния до фотонного детектора"<<endl;
		return 1;
	}
	if( ri1*opbeta < 1.0 ) {
		cerr<<"ri1="<<ri1<<" beta="<<opbeta
			<<": допороговая скорость в первом слое"<<endl;
		return 1;
	}
	if( optind >= argc )
		Usage(1);

	qefn=argv[optind];

//==========Расчет========//
	cout<<"________________________________________________\n"
		<<"Вычисление разрешения черенковского детектора.\n"
		<<"Параметры:\n"
		<<"  файл с квантовой эффективностью:\n   "<<qefn<<"\n"
		<<"  фактор эффективности ФЭУ:                     "<<efficiency<<"\n"
		<<"  размер пикселя ФЭУ:                           "<<pixelsize<<" мм\n"
		<<"  расстояние от радиатора до ФЭУ:               "<<D<<" мм\n"
		<<"  показатель преломления первого слоя на 400нм: "<<ri1<<"\n"
		<<"  число слоев радиатора:                        "<<nlayers<<"\n"
		<<"  полная толщина радиатора:                     "<<T<<" мм\n"
		<<"  длина рассеяния в аэрогеле на 400 нм:         "<<Lsc<<" мм\n"
		<<"  оптимизация для скорости:                     "<<opbeta<<"\n"
		<<"  минимизация "<<(minimize?"включена":"выключена")<<"\n";
	if( *outfn ) cout<<"  выходной root-файл: "<<outfn<<"\n";
	if( *macfn ) cout<<"  файл макроса для Geant4: "<<macfn<<"\n";
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
	cout<<"Квантовая эффективность определена от "<<wl1<<" до "<<wl2<<" нм, всего "
		<<phdeteff.Size()<<" значений."<<endl;
	cout<<"Длина волны максимальной чувствительности к ЧИ: "<<wl0<<" нм"<<endl;
	cout.precision(0);

	//Создаем многослойный аэрогелевый радиатор по идеальной модели без дисперсии для
	//длины волны максимальной чувствительности
	double t0=D-T; //proximity distance
	MLADescription mla0(t0,opbeta,wl0);
	mla0.SetScatteringLength(Lsc);

	mla0.MakeFixed(nlayers,D,ri1);
//	for(int l=0; l<nlayers; l++)
//		mla0.AddAlayer(ri1,T/nlayers);

	cout<<"Исходный аэрогелевый радиатор:"<<endl;
	mla0.Print("  ");

	MLAResult& res=mla0.Calculate(phdeteff,pixelsize);

	cout<<"  Среднее число фотоэлектронов: "<<res.npe<<endl;
	cout<<"  Средний радиус:               "<<res.radius<<endl;
	cout<<"  Ошибка радиуса на 1 фотон:    "<<res.sigma1<<endl;
	cout<<"  Ошибка радиуса на трек:       "<<res.sigma_t<<endl;

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

		cout<<"Первая попытка"<<endl;

		ROOT::Minuit2::FunctionMinimum minimum = combmin(1000,0.1);

		cout<<minimum;

		ROOT::Minuit2::MnUserParameterState state=minimum.UserState();

		if( !minimum.IsValid() ) {
			cout<<"Не удалось оптимизировать 1 :("<<endl;
			return 0;
		}

		cout<<"Освобождаем параметры"<<endl;
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
			cout<<"Не удалось оптимизировать 2 :("<<endl;
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

		cout<<"Аэрогелевый радиатор оптимизированный"<<endl;
		mla.Print("  ");

		res=mla.Calculate(phdeteff,pixelsize);
		cout<<"  Среднее число фотоэлектронов: "<<res.npe<<endl;
		cout<<"  Средний радиус:               "<<res.radius<<endl;
		cout<<"  Ошибка радиуса на 1 фотон:    "<<res.sigma1<<endl;
		cout<<"  Ошибка радиуса на трек:       "<<res.sigma_t<<endl;
	}

	TH1D hrad("hrad","Radius photoelectron distribution;radius, mm",Nr,res.rmin,res.rmax);
	for(int i=0; i<Nr; i++) hrad.SetBinContent(i+1,res.s[i]);

	if( *outfn ) {
		cout<<"Сохраняем распределение фотонов по радиусу в "<<outfn<<endl;
		TFile *outfile=new TFile(outfn,"RECREATE");
		if( outfile->IsOpen() )
			hrad.Write("hrad");
		else
			cerr<<outfn<<": не могу открыть файл на запись"<<endl;
		outfile->Close();
		delete outfile;
	}

	if( *macfn ) {
		cout<<"Сохраняем описание RICH в командный файл Geant4 "<<macfn<<endl;
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
			cerr<<macfn<<": не могу открыть файл на запись"<<endl;
		}
	}

	if( !batch ) {
		TCanvas *c=new TCanvas("c1","Распределение фотоэлектронов по радиусу");

		hrad.Draw("l");

		app->Run(kTRUE);
	}

	return 0;
}


