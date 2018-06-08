#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <libgen.h>
#include <string>
#include <set>

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
 
#include "MLADescription.h"
#include "Spectrum.h"

using namespace std;

static const char *optstring="qn:s:N:D:T:p:b:e:o:O:m::";
static const char *progname;

//Параметры программы
static string qefn;
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
static string mintype="Minuit2";
static string minalgo="Migrad";
static string outfn = "hrad_t.root";
static string macfn;

static MLADescription mla0, mlamin;
static Spectrum phdeteff;
static double penaltyCoef = 1.0;

static double get_max_sensitivity_wl()
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

static const set<string>& define_minimizer_options() 
{
    set<string> &defmintypes = *(new set<string>);
    defmintypes.insert("Minuit");
    defmintypes.insert("Minuit2");
    defmintypes.insert("GSLMultiMin");
    defmintypes.insert("GSLMultiFit");
    defmintypes.insert("GSLSimAn");
    defmintypes.insert("Linear");
    defmintypes.insert("Fumili");
    defmintypes.insert("Genetic");
    return defmintypes;
}

static double resolution(const double* x)
{
    int nlayers = mlamin.GetNlayers();
  	double T = mlamin.GetTotalThickness();
	double Teval = 0.;
    for(int l=0; l<nlayers; l++) {
	    if( l>0 ) mlamin.SetIndex(l,x[2*l-1]);
	    mlamin.SetThickness(l,x[2*l]);
	    Teval += x[2*l];
    }

	MLAResult& res=mlamin.Calculate(phdeteff, pixelsize);

	if( !res.valid ) return NAN;
	
	double penalty = Teval<=T?0.:penaltyCoef*(T-Teval);

	return (res.sigma_t+penalty)/mlamin.GetProximityDistance();
}

static void write_geant4_macfile(string macfn,MLADescription& mla)
{
    cout<<"Сохраняем описание RICH в командный файл Geant4 "<<macfn<<endl;
	ofstream macfile(macfn.c_str());
	if( macfile.is_open() ) {
		macfile<<"/Rich/pmt/qeDataFile "<<qefn<<"\n"
			<<"/Rich/pmt/detection "<<efficiency<<"\n"
			<<"/Rich/radiator/scatterLength "<<Lsc<<" mm\n"
			<<"/Rich/radiator/absDataFile none\n"
			<<"/Rich/mode manual\n"
			<<"/Rich/proximity "<<mla.GetProximityDistance()<<" mm\n";
		for(int l=0; l<nlayers; l++) {
			macfile<<"/Rich/radiator/addLayer "<<setprecision(5)<<mla.GetIndex(l)
				<<" "<<mla.GetThickness(l)<<" mm\n";
		}
		macfile<<"/Rich/update"<<endl;
	} else {
		cerr<<macfn<<": не могу открыть файл на запись"<<endl;
	}
}

void Usage(int status)
{
	cout<<"Usage: "<<progname<<" [OPTIONS] qefile\n"
		<<"Вычисляет разрешение детектора черенковских колец с многослойным фокусирующим "
		<<"аэрогелем. Оптимизирует разрешение. Параметр qefile задает путь к файлу с "
		<<"данными о квантовой эффективности фотонного детектора\n"
		<<" OPTIONS:\n"
		<<"   -q                Задачный режим без графики и интерпретатора\n"
		<<"   -e eff            Фактор к эффективности фотонного детектора 0<eff<=1 ("<<efficiency<<")\n"
		<<"   -p size           Размер квадратного пикселя фотонного детектора, мм ("<<pixelsize<<")\n"
		<<"   -D distance       Расстояние от начала радиатора до фотонного детектора, мм ("<<D<<")\n"
		<<"   -n ri1            Максимальный показатель преломления радиатора на 400 нм ("<<ri1<<")\n"
		<<"   -N nlayers        Число слоев аэрогеля ("<<nlayers<<")\n"
		<<"   -T thickness      Толщина радиатора, мм ("<<T<<")\n"
		<<"   -s Lsc            Длина рассеяния на 400 нм, мм ("<<Lsc<<")\n"
		<<"   -b beta           Оптимизировать радиатор для данной скорости ("<<opbeta<<")\n"
		<<"   -m [type[:algo]]  Оптимизировать радиатор по ошибке на трек с данным типом минимизатора и алгоритмом ("<<mintype<<":"<<minalgo<<")\n"
		<<"   -o filename       Выходной root файл для сохранения распределения по радиусу ("<<outfn<<")\n"
		<<"   -O filename       Сохранить описание детектора в командный файл Geant4 ("<<macfn<<")\n"
		<<endl;
	exit(status);
}

int main(int argc, char* argv[])
{
	progname=argv[0];

	if( argc==1 ) Usage(0);

    const set<string> defmintypes = define_minimizer_options();
    
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
			if( optarg ) {
			    string minarg = optarg;
			    size_t found = minarg.find_first_of(':');
			    if( found!=string::npos ) {
    			    mintype = minarg.substr(0,found-1);
    			    minalgo = minarg.substr(found+1);
			    } else {
    			    mintype = minarg;
    			    minalgo = "";
    			}
			    if( defmintypes.count(mintype)==0 ) {
			        cerr<<mintype<<": неизвестный тип минимизатора"<<endl;
			        return 1;
			    }
			}
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
		if( minimize )
    		cout << "     тип:             " << mintype << "\n"
    		     << "     алгоритм:        " << minalgo << "\n";
	if( !outfn.empty() ) cout<<"  выходной root-файл: "<<outfn<<"\n";
	if( !macfn.empty() ) cout<<"  файл макроса для Geant4: "<<macfn<<"\n";
	cout<<"________________________________________________"<<endl;

	TRint *app=0;
	if( !batch ) {
		int ac=2;
		char *av[2]={"rint","-l"};
		app=new TRint("rint",&ac,av);
	}

	phdeteff.ReadFile(qefn.c_str());
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

	cout<<"Исходный аэрогелевый радиатор:"<<endl;
	mla0.Print("  ");

	MLAResult& res=mla0.Calculate(phdeteff,pixelsize);

	cout<<"  Среднее число фотоэлектронов: "<<res.npe<<endl;
	cout<<"  Средний радиус:               "<<res.radius<<endl;
	cout<<"  Ошибка радиуса на 1 фотон:    "<<res.sigma1<<endl;
	cout<<"  Ошибка радиуса на трек:       "<<res.sigma_t<<endl;

	if( minimize ) {
		cout << "Оптимизация радиатора (тип: " << mintype << ", алгоритм: " << minalgo << ")" << endl;
        ROOT::Math::Minimizer* min = 
              ROOT::Math::Factory::CreateMinimizer(mintype, minalgo);
        min->SetMaxFunctionCalls(1000000);
        min->SetMaxIterations(100000);
        min->SetTolerance(0.001);
        min->SetPrintLevel(1);
     
        ROOT::Math::Functor fcn(&resolution,2*nlayers-1); 

		mlamin = mla0;
     
        for(int i=0; i<5; i++, penaltyCoef*=2) {
            min->SetFunction(fcn);
            // Set the free variables to be minimized
		    char name[20];
            min->SetLimitedVariable(0, "t1", mlamin.GetThickness(0), 0.01, 0.0, T);
		    for(int l=1; l<nlayers; l++) {
			    sprintf(name,"n%d",l+1);
                min->SetLimitedVariable(2*l-1, name, mlamin.GetIndex(l), 0.005, 1.0, 1.2);
			    sprintf(name,"t%d",l+1);
                min->SetLimitedVariable(2*l, name, mlamin.GetThickness(l), 0.5, 0.0, T);
		    }

		    cout<<"Попытка "<<i+1<<", penaltyCoef="<<penaltyCoef<<endl;
            min->Minimize(); 

            //min->PrintResults();

		    if( min->Status()!=0 ) {
			    cout<<">> Не удалось оптимизировать :("<<endl;
			    break;
		    }

            const double *xs = min->X();
		    for(int l=0; l<nlayers; l++) {
			    if( l>0 ) mlamin.SetIndex(l,xs[2*l-1]);
			    mlamin.SetThickness(l,xs[2*l]);
		    }
		    if( i<4 ) min->Clear();
		}
		mla0 = mlamin;

		cout<<"Аэрогелевый радиатор оптимизированный"<<endl;
		mla0.Print("  ");

		res=mla0.Calculate(phdeteff,pixelsize);
		cout<<"  Среднее число фотоэлектронов: "<<res.npe<<endl;
		cout<<"  Средний радиус:               "<<res.radius<<endl;
		cout<<"  Ошибка радиуса на 1 фотон:    "<<res.sigma1<<endl;
		cout<<"  Ошибка радиуса на трек:       "<<res.sigma_t<<endl;
	}

	TH1D hrad("hrad","Radius photoelectron distribution;radius, mm",Nr,res.rmin,res.rmax);
	for(int i=0; i<Nr; i++) hrad.SetBinContent(i+1,res.s[i]);

	if( !outfn.empty() ) {
		cout<<"Сохраняем распределение фотонов по радиусу в "<<outfn<<endl;
		TFile *outfile=new TFile(outfn.c_str(),"RECREATE");
		if( outfile->IsOpen() )
			hrad.Write("hrad");
		else
			cerr<<outfn<<": не могу открыть файл на запись"<<endl;
		outfile->Close();
		delete outfile;
	}

	if( !macfn.empty() )
		write_geant4_macfile(macfn, mla0);

	if( !batch ) {
		TCanvas *c=new TCanvas("c1","c1");

		hrad.Draw("l");

		app->Run(kTRUE);
	}

	return 0;
}


