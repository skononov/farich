const Double_t M[2]={139.57 ,493.6}; //MeV
const Double_t Mpi=M[0], Mk=M[1];
Double_t Lsc=40; //mm
Double_t A=1;
Double_t wl0=400.0; //wavelength, nm

Double_t WLmin=200, WLmax=700;

struct {
	Double_t t, n, p, d, g;
} InPar;

struct Result {
	Double_t Npe[2];
	Double_t Pth[2];
	Double_t Ang[2],  sAngAir[2], AngAir[2], dAngAir;
	Double_t sR[2][3], sRmr[3];
	Double_t Sep;
} r;

TSpline3 *sqe, *sqdis, *snpe[2];
TGraph *gqe, *gqdis, *gnpe;


Double_t func_qe(Double_t *x,Double_t *par)
{
	return sqe->Eval(*x);
}

Double_t func_index(Double_t *x,Double_t *par)
{
	return sqrt(1+(par[0]**2-1)*(sqdis->Eval(*x)**2-1)/(sqdis->Eval(wl0)**2-1));
}

Double_t func_ang_vs_wl(Double_t *x,Double_t *par)
{
	Double_t m=par[0];

	return acos(sqrt(1+(m/InPar.p)**2)/func_index(x,&par[1]));
}

Double_t func_angair_vs_wl(Double_t *x,Double_t *par)
{
	return asin(sin(func_ang_vs_wl(x,par))*func_index(x,&par[1]));
}

Double_t func_npe_vs_wl(Double_t *x,Double_t *par)
{
	Double_t rwl=(*x)/wl0, val, angle;

	angle=func_ang_vs_wl(x,par);

	val = 2*PI/137*A*InPar.g*Lsc/(wl0*1e-6)*sin(angle)**2*cos(angle)/wl0;
	val *= rwl**2*func_qe(x,0)/100*(1-exp(-InPar.t/rwl**4/Lsc/cos(angle)));

	return val;
}

Double_t func_npe_pi_vs_ang(Double_t *x,Double_t *par)
{
	return snpe[0]->Eval(*x);
}
Double_t func_npe_K_vs_ang(Double_t *x,Double_t *par)
{
	return snpe[1]->Eval(*x);
}

Int_t getres(Double_t t, Double_t pitch, Double_t n, Text_t* qefn, Double_t p, Double_t d, Double_t g=1.0, Bool_t out=false)
{
//	if( !gInterpreter->IsLoaded("dataread.C") ) gROOT->LoadMacro("dataread.C");
//	if( !gInterpreter->IsLoaded("getwidth.C") ) gROOT->LoadMacro("getwidth.C");

	InPar.t=t;
	InPar.n=n;
	InPar.p=p;
	InPar.d=d;
	InPar.g=g;

	for(int i=0; i<2; i++) {
		r.Npe[i]=0;
		r.Ang[i]=0;
		r.sAngAir[i]=Nan;
		r.AngAir[i]=0;
		r.sR[i][0]=Nan;
		r.sR[i][1]=Nan;
		r.sR[i][2]=Nan;
		r.sRmr[0]=Nan;
		r.sRmr[1]=Nan;
		r.sRmr[2]=Nan;
	}

	if( out ) cout<<"\n\n";

	TString qepath = "/home/skononov/work/rich/qe/";
	qepath += qefn;
    qepath += ".dat";

	gqe = new TGraph(qepath);
	WLmin = (gqe->GetX())[0];
	WLmax = (gqe->GetX())[gqe->GetN()-1];

	gqdis = new TGraph("/home/skononov/work/rich/quartzdis.dat");

	sqe = new TSpline3("QE",gqe);
	sqdis = new TSpline3("Quartz dispersion",gqdis);

	TF1 *Fqe = new TF1("Fqe",func_qe,WLmin,WLmax,0);
	TF1 *Find = new TF1("Find",func_index,WLmin,WLmax,1);
	TF1 *Fang = new TF1("Fang",func_ang_vs_wl,WLmin,WLmax,2);
	TF1 *Fangair = new TF1("Fangair",func_angair_vs_wl,WLmin,WLmax,2);
	TF1 *Fnpe_wl = new TF1("Fnpe_wl",func_npe_vs_wl,WLmin,WLmax,2);
    Fnpe_wl->SetTitle("Spectrum of detected photoelectrons");
	TF1 *Fnpe_ang[2];
	Find->SetParameter(0,n);
	Fang->SetParameters(Mpi,n);
	Fangair->SetParameters(Mpi,n);
	Fnpe_wl->SetParameters(Mpi,n);

	if (Find->Eval(WLmax)<sqrt(1+(Mpi/p)**2)) {//Pions below threshold in aerogel
		cout<<"Pions are below threshold. No separation.\n";
		return 1;
	}

	Bool_t KaboveThres = kTRUE;
	if (Find->Eval(WLmax)<sqrt(1+(Mk/p)**2)) //Kaons below threshold in aerogel
		KaboveThres=kFALSE;

	const Int_t np=31;
	Double_t ang[np], npe[np], wl[np], AngMin[2], AngMax[2];

	for(int i=0; i<2; i++)
	{
		if (i==1&&!KaboveThres) break;
		Fangair->SetParameter(0,M[i]);
		Fnpe_wl->SetParameter(0,M[i]);
		for(int ip=0; ip<np; ip++)
		{
			wl[ip] = WLmin+ip*(WLmax-WLmin)/(np-1);
			ang[np-ip-1] = 1e3*Fangair->Eval(wl[ip]); //mrad
			npe[np-ip-1] = TMath::Abs(Fnpe_wl->Eval(wl[ip])/Fangair->Derivative(wl[ip])/1e3);
		}
		if( gnpe ) delete gnpe;
		gnpe = new TGraph(np,ang,npe);
		snpe[i] = new TSpline3("Npe vs angle",gnpe);
		AngMin[i] = ang[TMath::LocMin(np,ang)];
		AngMax[i] = ang[TMath::LocMax(np,ang)];
	}
	Fnpe_ang[0] = new TF1("Fnpe_ang_pi",func_npe_pi_vs_ang,AngMin[0],AngMax[0],0);
    Fnpe_ang[0]->SetTitle("Angular density of photoelectrons");
	Fnpe_ang[1] = new TF1("Fnpe_ang_K",func_npe_K_vs_ang,AngMin[1],AngMax[1],0);
    Fnpe_ang[1]->SetTitle("Angular density of photoelectrons");

	Double_t rms_ch[2]={0,0};
	for(Int_t i=0; i<2; i++) {
		if (i==1&&!KaboveThres) break;
		Fnpe_wl->SetParameter(0,M[i]);
		r.Npe[i] = Fnpe_wl->Integral(WLmin,WLmax);
		r.Pth[i] = M[i]/sqrt(n**2-1);
		Fang->SetParameter(0,M[i]);
		Fangair->SetParameter(0,M[i]);
		Double_t WLmean=Fnpe_wl->Mean(WLmin,WLmax);
		r.Ang[i] = Fang->Eval(WLmean);
		r.AngAir[i] = Fangair->Eval(WLmean);
		r.sR[i][0] = t*tan(r.Ang[i]);
		r.sR[i][1] = pitch;
                rms_ch[i]=sqrt(Fnpe_ang[i]->Moment(2,AngMin[i],AngMax[i])-Fnpe_ang[i]->Mean(AngMin[i],AngMax[i])**2);
                Double_t cos2=cos(r.AngAir[i])**2;
		r.sR[i][2] = 1e-3*rms_ch[i]*sqrt(12)*d/cos2;
                r.sRmr[0] = 1e3*r.sR[0][0]/d*cos2/sqrt(12);
                r.sRmr[1] = 1e3*r.sR[0][1]/d*cos2/sqrt(12);
                r.sRmr[2] = rms_ch[0];

		Float_t Npes=r.Npe[i]>1?r.Npe[i]:1;
		r.sAngAir[i] = sqrt((r.sR[i][0]**2+r.sR[i][1]**2+r.sR[i][2]**2)/Npes/12)*cos(r.AngAir[i])**2/d;
	}
	r.dAngAir = r.AngAir[0]-r.AngAir[1];
	if (KaboveThres)
		r.Sep = r.dAngAir/sqrt(r.sAngAir[0]**2+r.sAngAir[1]**2);
	else
		r.Sep = -gausin(exp(-r.Npe[0]));

	if( out )
	{
		cout<<"--------------------Input data-----------------------------------\n";
		printf("RI=%1.3f, Pth(pi)=%4.1f MeV/c, Pth(K)=%4.1f MeV/c\n",n,r.Pth[0],r.Pth[1]);
		printf("Thickness=%2.1f mm, Lsc=%2.1f mm, A=%1.2f\n",t,Lsc,A);
		cout<<"-----------------------------------------------------------------\n";
		printf("Momentum %5.1fMeV/c\n",p);
		printf("Angle(pi)=%3.1f mrad, Angle(K)=%3.1f mrad, difference=%3.1f mrad\n",
			1e3*r.Ang[0], 1e3*r.Ang[1], 1e3*(r.Ang[0]-r.Ang[1]));
		printf("AngleAir(pi)=%3.1f mrad, AngleAir(K)=%3.1f mrad, difference=%3.1f mrad\n",
			1e3*r.AngAir[0], 1e3*r.AngAir[1], 1e3*r.dAngAir);
		printf("Level(pi)=%1.3f, Level(K)=%1.3f\n",1-(r.Pth[0]/p)**2,1-(r.Pth[1]/p)**2);
		cout<<"-----------------------------------------------------------------\n";
		cout<<"QE file name: "<<qepath<<endl;
		cout<<"Wavelength region: "<<WLmin<<"-"<<WLmax<<endl;
		printf("Proximity distance %3.0f mm, Photoelectron collection %1.2f\n",d,g);
		cout<<"-----------------------------------------------------------------\n";
        cout<<"************************RESULTS**********************************\n";
		cout<<"Npe(pi)="<<r.Npe[0]<<", Npe(K)="<<r.Npe[1]<<endl;
		cout<<"Npe_ang(pi)="<<Fnpe_ang[0]->Integral(AngMin[0],AngMax[0])<<", Npe_ang(K)="<<Fnpe_ang[1]->Integral(AngMin[1],AngMax[1])<<endl;
		printf("Ang.res.(pi)=%3.1f mrad, Ang.res.(K)=%3.1f mrad\n",1e3*r.sAngAir[0],1e3*r.sAngAir[1]);
		printf("Separation power %2.1f\n",r.Sep);
		printf("Contributions into the resolution, mm: Thickness  Pixel   Dispersion\n");
		printf("Pions   \t\t\t\t %2.1f       %2.1f      %2.1f\n", r.sR[0][0],r.sR[0][1],r.sR[0][2]);
		printf("Kaons   \t\t\t\t %2.1f       %2.1f      %2.1f\n", r.sR[1][0],r.sR[1][1],r.sR[1][2]);
		printf("Contributions to resolution, mrad: \t %2.1f       %2.1f      %2.1f\n",r.sRmr[0],r.sRmr[1],r.sRmr[2]);
		cout<<"*****************************************************************\n\n\n";
	}
/*	delete sqe;
	delete sqdis;
	delete snpe[0];
	delete snpe[1];
	delete gqe;
	delete gqdis;
	delete gnpe;
*/
	return 0;
}

