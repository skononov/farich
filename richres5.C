#include <iostream>
#include <iomanip>
#include "TSpline.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <Math/DistFunc.h>
#include "Math/SpecFunc.h"

//                        e      mu      pi      K       p
const Double_t Mtab[5]={0.511, 105.66, 139.57, 493.68, 938.27}; //MeV
Double_t Lsc=50; //mm
Double_t A=0.95;
Double_t wl0=400.0; //wavelength, nm

Double_t WLmin=200, WLmax=700;

const Double_t eulergamma=0.5772156649;

struct InputParameters {
    Double_t t, n, p, d, g;
} inPar;

struct Result {
    Double_t Npe[2];
    Double_t Pth[2];
    Double_t Ang[2],  AngAir[2], sAngAirSP[2], sAngAirT[2], dAngAir;
    Double_t Rad[2], sRadSP[2], sRadT[2];
    Double_t sRadSPc[2][3], sAngSPc[3];
    Double_t Sep;
} r;

TSpline3 *sqe, *snpe[2];
TGraph *gqe;


Double_t func_qe(Double_t *x,Double_t *par)
{
    return sqe->Eval(*x);
}

/*Refractive index as function of wavelength
    for aerogel of the nominal ref. ind. at 400nm.
  Scaled from LHCb data:
  T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
  Eur. Phys. J. C 52 (2007) 759-764
*/
Double_t func_aer_index(Double_t *x,Double_t *par)
{
    const Double_t LHCb_a0=0.05639, LHCb_wl0sqr = 83.22*83.22;
    Double_t wl = x[0];
    Double_t n400 = par[0];
    Double_t LHCb_RI2m1ref = LHCb_a0/(1-LHCb_wl0sqr/(400*400)); //(n**2-1) of LHCb aerogel at 400nm
    Double_t ri2m1_lhcb = LHCb_a0/(1-LHCb_wl0sqr/(wl*wl)); //(n**2-1) of LHCb aerogel at wl

    return sqrt( 1 + (n400*n400-1)/LHCb_RI2m1ref*ri2m1_lhcb );
}

// Original data: Ciddor 1996: n 0.23-1.690 µm, https://doi.org/10.1364/AO.35.001566
// by https://refractiveindex.com
Double_t func_air_index(Double_t *x,Double_t *par)
{
    const Double_t c[4] = {0.05792105, 238.0185, 0.00167917, 57.362};
    Double_t rwlum2 = 1./(1e-6*x[0]*x[0]);

    return 1 + c[0]/(c[1]-rwlum2) + c[2]/(c[3]-rwlum2);
}

Double_t func_ang_vs_wl(Double_t *x,Double_t *par)
{
    Double_t m=par[0];

    return acos(sqrt(1+pow((m/inPar.p),2))/func_aer_index(x,&par[1]));
}

Double_t func_anginair_vs_wl(Double_t *x,Double_t *par)
{
    return asin(sin(func_ang_vs_wl(x,par))*func_aer_index(x,&par[1])/func_air_index(x,NULL));
}

Double_t func_aer_npe_vs_wl(Double_t *x,Double_t *par)
{
    Double_t rwl=(*x)/wl0, val, angle;

    angle=func_ang_vs_wl(x,par);

    val = 2*M_PI/137*A*inPar.g*Lsc/(wl0*1e-6)*pow(sin(angle),2)*cos(angle)/wl0;
    val *= rwl*rwl*func_qe(x,0)/100*(1-exp(-inPar.t/pow(rwl,4)/Lsc/cos(angle)));

    return val;
}

Double_t func_air_npe_vs_wl(Double_t *x,Double_t *par)
{
    Double_t wl=*x, val, angle;

    angle=acos(sqrt(1+pow((par[0]/inPar.p),2))/func_air_index(x,NULL));

    val = 2*M_PI/137*inPar.g*inPar.d/(wl*1e-6)*pow(sin(angle),2)/wl*func_qe(x,0)/100;

    return val;
}

Double_t func_npe_vs_ang(Double_t *x,Double_t *par)
{
    int ipart=int(par[0]+0.5);
    if( !snpe[ipart] ) return 0;
    if( (*x)<snpe[ipart]->GetXmin() || (*x)>snpe[ipart]->GetXmax() ) return 0;
    return snpe[ipart]->Eval(*x);
}


Double_t separation_poisson(Double_t *x,Double_t *par)
{
    UInt_t nthr = int(*x+0.5);
    Double_t mu1 = par[0], mu2 = par[1];
    
    Double_t P1 = ROOT::Math::poisson_cdf(nthr,mu1), P2 = ROOT::Math::poisson_cdf(nthr,mu2);
    Double_t Q1 = P1==1.?0.:ROOT::Math::normal_quantile(P1,1), Q2 = P2==1.?0.:ROOT::Math::normal_quantile(P2,1);
    
    return Q2-Q1;
}

/** Calculates angular and radius resolution and number of photoelectrons for an aerogel proximity focusing RICH
    Parameters:
        t       - thickness of aerogel, mm
        pixel   - PD pixel size, mm
        n       - refractive index of aerogel
        qefn    - QE data file name without .dat extension. It should be located in ../qe/ directory.
        p       - particle momentum, MeV/c
        d       - air gap size between aerogel and photon detector, mm
        g       - additional efficiency coefficient
        part1,part2 - particle name: e, mu, pi, K, p
        out     - if kTRUE, print results
*/
Int_t getres(Double_t t, Double_t pixel, Double_t n, const Char_t* qefn, Double_t p,
             Double_t d, Double_t g=1.0, TString part1="pi", TString part2="K", Bool_t out=kTRUE)
{
    inPar.t=t;
    inPar.n=n;
    inPar.p=p;
    inPar.d=d;
    inPar.g=g;

    r.dAngAir=0;
    for(int i=0; i<2; i++) {
        r.Npe[i]=0;
        r.Ang[i]=0;
        r.AngAir[i]=0;
        r.sAngAirSP[i]=-1;
        r.sAngAirT[i]=-1;
        r.sRadSP[i]=-1;
        r.sRadT[i]=-1;
        r.sRadSPc[i][0]=-1;
        r.sRadSPc[i][1]=-1;
        r.sRadSPc[i][2]=-1;
        r.sAngSPc[0]=-1;
        r.sAngSPc[1]=-1;
        r.sAngSPc[2]=-1;
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
            cerr << "Unknown particle ID: " << spart << endl;
            return 1;
        }
    }

    if( out ) cout << "\n\n";

    if( M[1]<M[0] ) {
        Double_t temp;
        temp=M[0]; M[0]=M[1]; M[1]=temp;
    }

    gqe = new TGraph(qefn);
    if( gqe->GetN()==0 ) return 2;
    WLmin = (gqe->GetX())[0];
    WLmax = (gqe->GetX())[gqe->GetN()-1];

    sqe = new TSpline3("QE",gqe);

    TF1 *Fqe = new TF1("Fqe",func_qe,WLmin,WLmax,0);
    TF1 *Find = new TF1("Find",func_aer_index,WLmin,WLmax,1);
    TF1 *Fairind = new TF1("Fairind",func_air_index,WLmin,WLmax,0);
    TF1 *Fang = new TF1("Fang",func_ang_vs_wl,WLmin,WLmax,2);
    TF1 *Fangair = new TF1("Fangair",func_anginair_vs_wl,WLmin,WLmax,2);
    TF1 *Fnpe_wl = new TF1("Fnpe_wl",func_aer_npe_vs_wl,WLmin,WLmax,2);
    TF1 *Fairnpe_wl = new TF1("Fairnpe_wl",func_air_npe_vs_wl,WLmin,WLmax,1);
    TF1 *Fsep_npe = new TF1("Fsep_npe",separation_poisson,0,100,2);
    Fnpe_wl->SetTitle("Spectrum of detected unscattered photoelectrons from aerogel;#lambda, nm;dN_{pe}}/d#lambda");
    Fairnpe_wl->SetTitle("Spectrum of detected photoelectrons from air;#lambda, nm;dN_{pe}/d#lambda");
    Find->SetParameter(0,n);
    Fang->SetParameters(M[0],n);
    Fangair->SetParameters(M[0],n);
    Fnpe_wl->SetParameters(M[0],n);
    Fairnpe_wl->SetParameter(0,M[0]);

    if (Find->Eval(WLmax)<sqrt(1+pow((M[0]/p),2))) {//Light particle below threshold in aerogel
        cout << "Light particle is below Cherenkov threshold. No separation.\n";
        return 3;
    }

    Bool_t p2aboveThres = kTRUE;
    if (Find->Eval(WLmax)<sqrt(1+pow((M[1]/p),2))) //Heavy particle below threshold in aerogel
        p2aboveThres=kFALSE;

    const Int_t np=31;
    Double_t ang[np], npe[np], AngMin[2] = {1e3*M_PI, 1e3*M_PI}, AngMax[2] = {0, 0}, airnpe[2] = {0, 0};

    for(int i=0; i<2; i++)
    {
        if (i==1 && !p2aboveThres) break;
        Fangair->SetParameter(0,M[i]);
        Fnpe_wl->SetParameter(0,M[i]);
        if (Fairind->Eval(WLmax) > sqrt(1+pow((M[i]/p),2)))
            airnpe[i] = Fairnpe_wl->Integral(WLmin, WLmax);
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
    TF1* Fnpe_ang = new TF1("Fnpe_ang",func_npe_vs_ang,TMath::MinElement(2,AngMin),TMath::MaxElement(2,AngMax),1);
    Fnpe_ang->SetTitle("Angular density of photoelectrons;angle, mrad;dN_{pe}/d#theta");

    Double_t rms_ch[2]={0,0};
    for(Int_t i=0; i<2; i++) {
        if (i==1 && !p2aboveThres) break;
        Fnpe_wl->SetParameter(0,M[i]);
        r.Npe[i] = Fnpe_wl->Integral(WLmin,WLmax);
        if( r.Npe[i]<1e-6 ) {
            cerr << "Too little photoelectrons for particle " << (i==0?part1:part2) << ". Close to threshold?" << endl;
            return 2;
        }
        r.Pth[i] = M[i]/sqrt(n*n-1);
        Fang->SetParameter(0,M[i]);
        Fangair->SetParameter(0,M[i]);
        Double_t WLmean=Fnpe_wl->Mean(WLmin,WLmax);
        r.Ang[i] = Fang->Eval(WLmean);
        r.AngAir[i] = Fangair->Eval(WLmean);
        r.sRadSPc[i][0] = t*tan(r.Ang[i])/sqrt(12);
        r.sRadSPc[i][1] = pixel/sqrt(12);
        Fnpe_ang->SetParameter(0,i);
        rms_ch[i] = sqrt(Fnpe_ang->Variance(Fnpe_ang->GetXmin(),Fnpe_ang->GetXmax()));
        Double_t cos2 = cos(r.AngAir[i])*cos(r.AngAir[i]);
        r.sRadSPc[i][2] = 1e-3*rms_ch[i]*(d+t/2)/cos2;

        Double_t factor = (ROOT::Math::expint(r.Npe[i])-eulergamma-log(r.Npe[i]))/(exp(r.Npe[i])-1);

        r.sRadSP[i] = sqrt(pow(r.sRadSPc[i][0],2)+pow(r.sRadSPc[i][1],2)+pow(r.sRadSPc[i][2],2));
        r.sRadT[i] = r.sRadSP[i]*sqrt(factor);
        r.sAngAirSP[i] = r.sRadSP[i]*cos2/d;
        r.sAngAirT[i] = r.sAngAirSP[i]*sqrt(factor);
        
        if( i==0 ) {
            r.sAngSPc[0] = 1e3*r.sRadSPc[0][0]/(d+t/2)*cos2;
            r.sAngSPc[1] = 1e3*r.sRadSPc[0][1]/(d+t/2)*cos2;
            r.sAngSPc[2] = rms_ch[0];
        }
    }
    r.dAngAir = r.AngAir[0]-r.AngAir[1];
    Double_t sep_res = p2aboveThres?2*r.dAngAir/(r.sAngAirT[0]+r.sAngAirT[1]):0.0;

    Fsep_npe->SetParameters(r.Npe[0], r.Npe[1]);
    Fsep_npe->SetRange(0, 2*r.Npe[0]);
    Double_t sep_npe = Fsep_npe->GetMaximum();

    Bool_t isSepByNpe = kFALSE;
    if (sep_npe > sep_res) {
        r.Sep = sep_npe;
        isSepByNpe = kTRUE;
    } else
        r.Sep = sep_res;

    if (out)
    {
        //		cout.setf(ios_base::scientific);
        cout << "--------------------------Input data------------------------------------\n"
             << setprecision(4) << "RI=" << n << ", Pth(" << part1 << ")=" << r.Pth[0] << " MeV/c, Pth("
             << part2 << ")=" << r.Pth[1] << " MeV/c\n"
             << setprecision(3) << "Thickness=" << t << " mm, Lsc=" << Lsc << " mm, A=" << A << "\n"
             << "----\n"
             << setprecision(4) << "P=" << p << " MeV/c\n"
             << "Angle(" << part1 << ")=" << 1e3*r.Ang[0] << " mrad, Angle(" << part2 << ")="
             << 1e3*r.Ang[1] << " mrad, difference=" << 1e3*(r.Ang[0]-r.Ang[1]) << " mrad\n"
             << "AngleAir(" << part1 << ")=" << 1e3*r.AngAir[0] << " mrad, AngleAir(" << part2 << ")="
             << 1e3*r.AngAir[1] << " mrad, difference=" << 1e3*r.dAngAir << " mrad\n"
             << "Intensity(" << part1 << ")=" << 1-pow(r.Pth[0]/p,2) << ", Intensity(" << part2 << ")="
             << 1-pow(r.Pth[1]/p,2) << "\n"
             << "----\n"
             << "QE file name: " << qefn << "\n"
             << "Wavelength region: " << WLmin << "-" << WLmax << "\n"
             << setprecision(3) << "Proximity distance " << d << " mm, Photoelectron collection " << g << "\n"
             << "------------------------------------------------------------------------\n" << endl;
        cout << "************************RESULTS*****************************************\n"
             << setprecision(3) << "Npe.wl(" << part1 << ")=" << r.Npe[0];
        if (p2aboveThres)
            cout << ", Npe.wl(" << part2 << ")=" << r.Npe[1];
        else
            cout << ", " << part2 << " is below threshold";
        Fnpe_ang->SetParameter(0,0);
        cout << "\nNpe.ang(" << part1 << ")=" << Fnpe_ang->Integral(AngMin[0],AngMax[0]);
        if (p2aboveThres) {
            Fnpe_ang->SetParameter(0,1);
            cout << ", Npe.ang(" << part2 << ")=" << Fnpe_ang->Integral(AngMin[1],AngMax[1]);
        } else
            cout << ", " << part2 << " is below threshold";
        if (airnpe[0]==0)
            cout << "\nBoth particles are below threshold in air";
        else {
            cout << "\nNpe.in_air(" << part1 << ")=" << airnpe[0];
            if (airnpe[1]>0)
                cout << ", Npe.in_air(" << part2 << ")=" << airnpe[1];
        }
        
        cout << "\n" << setprecision(3) << setw(5) << "Ang.res.per_track(" << part1 << ")=" << 1e3*r.sAngAirT[0] << " mrad, Ang.res.per_photon(" << part1 << ")=" << 1e3*r.sAngAirSP[0] << " mrad"
             << "\nRad.res.per_track(" << part1 << ")=" << r.sRadT[0] << " mm, Rad.res.per_photon(" << part1 << ")=" << r.sRadSP[0] << " mm";
        if (p2aboveThres) {
             cout << "\nAng.res.per_track(" << part2 << ")=" << 1e3*r.sAngAirT[1] << " mrad, Ang.res.per_photon(" << part2 << ")=" << 1e3*r.sAngAirSP[1] << " mrad"
                  << "\nRad.res.per_track(" << part2 << ")=" << r.sRadT[1] << " mm, Rad.res.per_photon(" << part2 << ")=" << r.sRadSP[1] << " mm";
        }
        cout << "\nSeparation power" << (isSepByNpe?"(by Npe!): ":": ") << setprecision(2) << r.Sep << "\n"
             << "Contributions into the resolution, mm:  Thickness   Pixel   Dispersion\n"
             << setw(2) << part1 << "\t\t\t\t\t" << setprecision(3)
             << setw(6) << r.sRadSPc[0][0] << "     " << setw(6) << r.sRadSPc[0][1] << "    " << setw(6) << r.sRadSPc[0][2] << "\n";
        if (p2aboveThres) {
            cout << setw(2) << part2 << "\t\t\t\t\t" << setprecision(3)
                 << setw(6) << r.sRadSPc[1][0] << "     " << setw(6) << r.sRadSPc[1][1] << "    " << setw(6) << r.sRadSPc[1][2] << "\n";
        } else {
            cout << setw(2) << part2 << "\t\t\t\t\t\tUnder threshold"  << "\n";
        }
        cout << "Contributions to resolution, mrad:\t" << setprecision(3)
             << setw(6) << r.sAngSPc[0] << "     " << setw(6) << r.sAngSPc[1] << "    " << setw(6) << r.sAngSPc[2] << "\n"
             << "************************************************************************\n" << endl;
        cout.precision(0);
        cout.width(0);
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

