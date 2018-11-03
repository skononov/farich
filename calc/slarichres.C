#include <iostream>
#include <iomanip>
//                        e      mu      pi      K       p
const Double_t Mtab[5]={0.511, 105.66, 139.57, 493.68, 938.27}; //MeV
Double_t Lsc=50; //mm
Double_t A=0.95;
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

TSpline3 *sqe, *snpe[2];
TGraph *gqe;


Double_t func_qe(Double_t *x,Double_t *par)
{
    return sqe->Eval(*x);
}

//Parameters of the one-pole Sellmeier formula fit to Novosibirsk aerogel n=1.03 data.
//T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
//Eur. Phys. J. C 52 (2007) 759-764
const Double_t LHCb_a0=0.05639, LHCb_wl0sqr = 83.22*83.22;
const Double_t LHCb_RI2m1ref = LHCb_a0/(1-LHCb_wl0sqr/(400*400));

Double_t func_index_lhcb(Double_t *x,Double_t *par)
{
        Double_t wl = *x; //wavelength at which calculate refractive index [nm]
        Double_t ri0 = *par; //given refractive index at wl0
        Double_t ri2m1_lhcb = LHCb_a0/(1-LHCb_wl0sqr/(wl*wl)); //refractive index of LHCb aerogel at wl
        return sqrt( 1 + (ri0*ri0-1)/LHCb_RI2m1ref*ri2m1_lhcb );
}

Double_t func_ang_vs_wl(Double_t *x,Double_t *par)
{
    Double_t m=par[0];

    return acos(sqrt(1+(m/InPar.p)**2)/func_index_lhcb(x,&par[1]));
}

Double_t func_angair_vs_wl(Double_t *x,Double_t *par)
{
    return asin(sin(func_ang_vs_wl(x,par))*func_index_lhcb(x,&par[1]));
}

Double_t func_npe_vs_wl(Double_t *x,Double_t *par)
{
    Double_t rwl=(*x)/wl0, val, angle;

    angle=func_ang_vs_wl(x,par);

    val = 2*PI/137*A*InPar.g*Lsc/(wl0*1e-6)*sin(angle)**2*cos(angle)/wl0;
    val *= rwl**2*func_qe(x,0)/100*(1-exp(-InPar.t/rwl**4/Lsc/cos(angle)));

    return val;
}

Double_t func_npe_vs_ang(Double_t *x,Double_t *par)
{
    int ipart=int(par[0]+0.5);
    if( (*x)<snpe[ipart]->GetXmin() || (*x)>snpe[ipart]->GetXmax() ) return 0;
    return snpe[ipart]->Eval(*x);
}

/* Calculate number of photoelectrons, radius&angle resolutions, 2-particle separation for
   single layer aerogel RICH.
   Input parameters:
    t           - thickness of radiator, mm
    pxsize      - PD pixel size, mm
    n           - refractive index of aerogel at 400nm
    qefn        - QE data file name (may be a file name in "../qe" directory without ".dat" extension
    p           - particle momentum in MeV/c
    d           - proximity distance (air gap size), mm
    g           - additional factor in efficiency 0..1
    part1,part2 - particle names to evaluate separation for (e,mu,pi,k,p), case insensitive
    out         - print out calculated values, otherwise - just store in 'struct Result r'
 */

Int_t getres(Double_t t, Double_t pxsize, Double_t n, const Char_t* qefn, Double_t p,
             Double_t d, Double_t g=1.0, TString part1="pi", TString part2="K", Bool_t out=kTRUE)
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

    Double_t M[2];
    for(int i=0; i<2; i++) {
        TString spart=(i==0?part1:part2);
        spart.ToLower();
        if( spart=="e" )
            M[i]=Mtab[0];
        else if( spart=="mu" )
            M[i]=Mtab[1];
        else if( spart=="pi" )
            M[i]=Mtab[2];
        else if( spart=="k" )
            M[i]=Mtab[3];
        else if( spart=="p" )
            M[i]=Mtab[4];
        else {
            cerr<<"Unknown particle ID: "<<spart<<endl;
            return 1;
        }
    }

    if( out ) cout<<"\n\n";

    if( M[1]<M[0] ) {
        Double_t temp;
        temp=M[0]; M[0]=M[1]; M[1]=temp;
    }

    TString qepath;
    if( gSystem->AccessPathName(qefn)==0 ) //file with QE data is found
        qepath = qefn;
    else {
        qepath = "data/";
        qepath += qefn;
        qepath += ".dat";
    }

    gqe = new TGraph(qepath);
    if( gqe->GetN()==0 ) return 2;
    WLmin = (gqe->GetX())[0];
    WLmax = (gqe->GetX())[gqe->GetN()-1];

    sqe = new TSpline3("QE",gqe);

    TF1 *Fqe = new TF1("Fqe",func_qe,WLmin,WLmax,0);
    TF1 *Find = new TF1("Find",func_index_lhcb,WLmin,WLmax,1);
    TF1 *Fang = new TF1("Fang",func_ang_vs_wl,WLmin,WLmax,2);
    TF1 *Fangair = new TF1("Fangair",func_angair_vs_wl,WLmin,WLmax,2);
    TF1 *Fnpe_wl = new TF1("Fnpe_wl",func_npe_vs_wl,WLmin,WLmax,2);
    Fnpe_wl->SetTitle("Spectrum of detected photoelectrons");
    Find->SetParameter(0,n);
    Fang->SetParameters(M[0],n);
    Fangair->SetParameters(M[0],n);
    Fnpe_wl->SetParameters(M[0],n);

    if (Find->Eval(WLmax)<sqrt(1+(M[0]/p)**2)) {//Light particle below threshold in aerogel
        cout<<"Light particle is below Cherenkov threshold. No separation.\n";
        return 3;
    }

    Bool_t KaboveThres = kTRUE;
    if (Find->Eval(WLmax)<sqrt(1+(M[1]/p)**2)) //Heavy particle below threshold in aerogel
        KaboveThres=kFALSE;

    const Int_t np=31;
    Double_t ang[np], npe[np], AngMin[2], AngMax[2];

    for(int i=0; i<2; i++)
    {
        if (i==1&&!KaboveThres) break;
        Fangair->SetParameter(0,M[i]);
        Fnpe_wl->SetParameter(0,M[i]);
        for(int ip=0; ip<np; ip++)
        {
            Double_t wl = WLmax-ip*(WLmax-WLmin)/(np-1);
            ang[ip] = 1e3*Fangair->Eval(wl); //mrad
            npe[ip] = TMath::Abs(Fnpe_wl->Eval(wl)/Fangair->Derivative(wl)/1e3);
        }
        snpe[i] = new TSpline3("Npe vs angle",ang,npe,np);
        AngMin[i] = ang[0];
        AngMax[i] = ang[np-1];
    }
    TF1* Fnpe_ang = new TF1("Fnpe_ang",func_npe_vs_ang,AngMin[1],AngMax[0],1);
    Fnpe_ang->SetTitle("Angular density of photoelectrons");

    Double_t rms_ch[2]={0,0};
    for(Int_t i=0; i<2; i++) {
        if (i==1&&!KaboveThres) break;
        Fnpe_wl->SetParameter(0,M[i]);
        r.Npe[i] = Fnpe_wl->Integral(WLmin,WLmax);
        if( r.Npe[i]<1e-6 ) {
            cerr<<"Too little photoelectrons for particle "<<(i==0?part1:part2)<<". Close to threshold?"<<endl;
            return 2;
        }
        r.Pth[i] = M[i]/sqrt(n**2-1);
        Fang->SetParameter(0,M[i]);
        Fangair->SetParameter(0,M[i]);
        Double_t WLmean=Fnpe_wl->Mean(WLmin,WLmax);
        r.Ang[i] = Fang->Eval(WLmean);
        r.AngAir[i] = Fangair->Eval(WLmean);
        r.sR[i][0] = t*tan(r.Ang[i])/sqrt(12);
        r.sR[i][1] = pxsize/sqrt(12);
        Fnpe_ang->SetParameter(0,i);
        rms_ch[i] = sqrt(Fnpe_ang->Variance(AngMin[i],AngMax[i]));
        Double_t cos2 = cos(r.AngAir[i])*cos(r.AngAir[i]);
        r.sR[i][2] = 1e-3*rms_ch[i]*(d+t/2)/cos2;
        Float_t Npes=r.Npe[i]>1?r.Npe[i]:1;
        r.sAngAir[i] = sqrt((r.sR[i][0]**2+r.sR[i][1]**2+r.sR[i][2]**2)/Npes)*cos2/d;
        if( i==0 ) {
            r.sRmr[0] = 1e3*r.sR[0][0]/(d+t/2)*cos2;
            r.sRmr[1] = 1e3*r.sR[0][1]/(d+t/2)*cos2;
            r.sRmr[2] = rms_ch[0];
        }
    }
    r.dAngAir = r.AngAir[0]-r.AngAir[1];
    if (KaboveThres)
        r.Sep = 2*r.dAngAir/(r.sAngAir[0]+r.sAngAir[1]);
    else
        r.Sep = -gausin(exp(-r.Npe[0]));

    if( out )
    {
        //		cout.setf(ios_base::scientific);
        cout<<"--------------------------Input data------------------------------------\n"
            <<setprecision(4)<<"RI="<<n<<", Pth("<<part1<<")="<<r.Pth[0]<<" MeV/c, Pth("
            <<part2<<")="<<r.Pth[1]<<" MeV/c\n"
            <<setprecision(3)<<"Thickness="<<t<<" mm, Lsc="<<Lsc<<" mm, A="<<A<<"\n"
            <<"----\n"
            <<setprecision(4)<<"P="<<p<<" MeV/c\n"
            <<"Angle("<<part1<<")="<<1e3*r.Ang[0]<<" mrad, Angle("<<part2<<")="
            <<1e3*r.Ang[1]<<" mrad, difference="<<1e3*(r.Ang[0]-r.Ang[1])<<" mrad\n"
            <<"AngleAir("<<part1<<")="<<1e3*r.AngAir[0]<<" mrad, AngleAir("<<part2<<")="
            <<1e3*r.AngAir[1]<<" mrad, difference="<<1e3*r.dAngAir<<" mrad\n"
            <<"Intensity("<<part1<<")="<<1-(r.Pth[0]/p)**2<<", Intensity("<<part2<<")="
            <<1-(r.Pth[1]/p)**2<<"\n"
            <<"----\n"
            <<"QE file name: "<<qepath<<"\n"
            <<"Wavelength region: "<<WLmin<<"-"<<WLmax<<"\n"
            <<setprecision(3)<<"Proximity distance "<<d<<" mm, Photoelectron collection "<<g<<"\n"
            <<"------------------------------------------------------------------------\n"<<endl;
        cout<<"************************RESULTS*****************************************\n"
            <<setprecision(3)<<"Npe("<<part1<<")="<<r.Npe[0]<<", Npe("<<part2<<")="<<r.Npe[1]<<"\n";
        Fnpe_ang->SetParameter(0,0);
        cout<<"Npe_ang("<<part1<<")="<<Fnpe_ang->Integral(AngMin[0],AngMax[0]);
        Fnpe_ang->SetParameter(0,1);
        cout<<", Npe_ang("<<part2<<")="<<Fnpe_ang->Integral(AngMin[1],AngMax[1])<<"\n";
        cout<<setprecision(3)<<setw(5)<<"Ang.res.("<<part1<<")="<<1e3*r.sAngAir[0]
            <<" mrad, Ang.res.("<<part2<<")="<<1e3*r.sAngAir[1]<<" mrad\n"
            <<"Separation power "<<setprecision(2)<<r.Sep<<"\n"
            <<"Contributions into the resolution, mm:  Thickness  Pixel   Dispersion\n"
            <<setw(2)<<part1<<"\t\t\t\t\t"<<setprecision(3)
            <<setw(6)<<r.sR[0][0]<<"     "<<setw(6)<<r.sR[0][1]<<"    "<<setw(6)<<r.sR[0][2]<<"\n"
            <<setw(2)<<part2<<"\t\t\t\t\t"<<setprecision(3)
            <<setw(6)<<r.sR[1][0]<<"     "<<setw(6)<<r.sR[1][1]<<"    "<<setw(6)<<r.sR[1][2]<<"\n"
            <<"Contributions to resolution, mrad:\t"<<setprecision(3)
            <<setw(6)<<r.sRmr[0]<<"     "<<setw(6)<<r.sRmr[1]<<"    "<<setw(6)<<r.sRmr[2]<<"\n"
            <<"************************************************************************\n"<<endl;
        cout.precision(0);
        cout.width(0);
    }
    return 0;
}

