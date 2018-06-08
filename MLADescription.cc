#include "MLADescription.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

#include "Spectrum.h"

using namespace std;

double MLADescription::tolerance = 1e-6;

MLADescription::MLADescription(double d, double b, double wl) :
    nlayers(0),
    t0(d),
    beta(b),
    wavelength(wl),
    scatteringLength(0),
    totalThickness(0),
    result()
{}

MLADescription::MLADescription(const MLADescription& mla) :
nlayers(mla.nlayers),
    t0(mla.t0),
    beta(mla.beta),
    wavelength(mla.wavelength),
    scatteringLength(mla.scatteringLength),
    vn(mla.vn),
    vt(mla.vt),
    totalThickness(mla.totalThickness),
    result(mla.result)
{}

MLADescription::~MLADescription()
{}

void MLADescription::clear()
{
    rn.clear();
    rt.clear();
    vn.clear();
    vt.clear();
    nlayers=0;
    totalThickness=0;
    result.valid=false;
}

void MLADescription::GoToAbs()
{
    nlayers=rn.size()-1;
    vn.resize(nlayers);
    vt.resize(nlayers);
    totalThickness=0;
    for(int i=0; i<nlayers; i++) {
        vn[i]=AerogelRefIndex(rn[i+1]/beta,wavelength,400.);
        vt[i]=rt[i+1]*t0;
        totalThickness+=vt[i];
    }
}

int MLADescription::MakeAlayer()
{
    int num = rn.size();

    rn.push_back(rn.back()); //zero approximation

    int nit=0;
    double S, disc;
    do { //iterate
        S = 0;
        for (int i=0; i<num; i++) {
            S += rt[i]*Skl(num,i);
        }
        disc = T10 - S*Tkl(num,num);
        S = T10/S;
        rn[num] = sqrt(1+S*S);
        nit++;
        if (nit>500) {
            rn.pop_back(); //reset this calculation
            return 1; //too many loops
        }
    } while ( fabs(disc) > T10*tolerance );

    rt.push_back(rt[1]*Tkl(1,1)/Tkl(num,num)); //derive thickness of the layer
    return 0;
}

/* Old aerogel RI calculation based on quartz refractive index
double MLADescription::AerogelRefIndex(double n, double wl1, double wl2)
{
    static const double c[6] = {0.6961663, 0.0046791, 0.4079426, 0.0135121, 0.8974794, 97.9340025};

    double x=1e-6*wl1*wl1; //wavelength^2 [um^2]
    double nq1=sqrt( 1 + x*(c[0]/(x-c[1]) + c[2]/(x-c[3]) + c[4]/(x-c[5])) );
    x=1e-6*wl2*wl2;
    double nq2=sqrt( 1 + x*(c[0]/(x-c[1]) + c[2]/(x-c[3]) + c[4]/(x-c[5])) );

    return sqrt( 1 + (nq2*nq2-1)*(n*n-1)/(nq1*nq1-1) );
}
*/

double MLADescription::AerogelRefIndex(double n, double wl1, double wl2)
{
//Parameters of the one-pole Sellmeier formula fit to Novosibirsk aerogel n=1.03 data.
//T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
//Eur. Phys. J. C 52 (2007) 759-764
    const double LHCb_a0=0.05639, LHCb_wl0sqr = 83.22*83.22;
    double LHCb_RI2m1ref = LHCb_a0/(1-LHCb_wl0sqr/(wl1*wl1)); //(n**2-1) of LHCb aerogel at wl1
    double ri2m1_lhcb = LHCb_a0/(1-LHCb_wl0sqr/(wl2*wl2)); //(n**2-1) of LHCb aerogel at wl2

    return sqrt( 1 + (n*n-1)/LHCb_RI2m1ref*ri2m1_lhcb );
}

void MLADescription::AddAlayer(double ri,double t)
{
    nlayers++;
    vn.push_back(ri);
    vt.push_back(t);
    totalThickness+=t;
}

int MLADescription::MakeLayers(int N, double n1, double t1)
{ //fixed number of layers. Input parameters: beta, t0, n1 (at 400nm), t1, N.
    clear();

    rn.reserve(N+1);
    rt.reserve(N+1);

    if( N<=0 ) return 0;

    if( wavelength!=400. ) n1=AerogelRefIndex(n1,400.,wavelength);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1*beta);
    rt.push_back(t1/t0);

    T10 = Tkl(1,0);

    int i=1;

    while (i < N) { //loop on layer number
        //derive index of the i-th layer
        if (MakeAlayer()!=0) break;
        i++;
    }

    GoToAbs();

    return i;
}

int MLADescription::MakeGabarit(double G, double n1, double rt1)
{ //fixed total thickness of radiator. Input parameters: beta, t0, G, n1 (at 400nm), rt1.

    double t1 = rt1*t0; //thickness of the first layer

    double Tmax=G-t0; //maximum thickness of radiator (may be unachievable)

    clear();

    if( wavelength!=400. ) n1=AerogelRefIndex(n1,400.,wavelength);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1*beta);
    rt.push_back(rt1);

    T10 = Tkl(1,0);

    double T=0, delta=t1;
    int i=0;

    do { //loop until required thickness gained
        T+=delta;
        if (MakeAlayer()!=0) return i; //something wrong with calculation
        delta=t0*rt.back(); //thickness of a latest layer evaluated
        i++;
    } while ( T+delta < Tmax );

    rt.pop_back();  //remove the last layer
    rn.pop_back();  //remove the last layer
    i--;

    GoToAbs();

    return i; //number of layers to fit given Tmax
}

int MLADescription::MakeFixed(int N, double G, double n1)
{ //fixed number of layers and total thickness of radiator.
    //Input parameters: beta, t0, N, G, n1 (at 400nm)
    clear();

    if( N<=0 ) return 0;

    if( N==1 ) { //single layer case
        nlayers=1;
        double t1=G-t0;
        vn.push_back(n1);
        vt.push_back(t1);
        totalThickness=t1;
        return 1;
    }

    if( wavelength!=400. ) n1=AerogelRefIndex(n1,400.,wavelength);

    rn.reserve(N+1);
    rt.reserve(N+1);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1*beta);
    T10 = Tkl(1,0);

    double Tgoal=G-t0; //total thickness of radiator to attain
    double T=0;

    double t1=Tgoal/N; //estimation of the 1-st layer thickness

    int nit=0;

    while ( fabs(T-Tgoal)>tolerance*Tgoal ) { //loop until required accuracy of thickness gained
        rn.resize(2);
        rt.resize(1);

        rt.push_back(t1/t0);
        T=t1;
        int i=1;
        for ( ; i<N; i++) {
            //derive the i-th layer
            if (MakeAlayer()!=0) break;
            T+=t0*rt[i+1];
        }
        if (i<N) break; //layers calculation failed

        t1+=(Tgoal-T)*t1/T;

        nit++;
        if (nit>500) return 0; //iterations diverging
    }

    GoToAbs();

    return nlayers;
}

void MLADescription::Print(const char* pfx) const
{
    if( !nlayers ) {
        cout<<"Radiator not defined yet!"<<endl;
        return;
    }

    cout.setf(ios::fixed|ios::left);
    cout.precision(2);
    cout << pfx << "Number of layers: " << nlayers << "\n"
        << pfx << "Proximity distance:  " << t0 << "\n"
        << pfx << "Total thickness:  " << totalThickness << "\n"
        << pfx << "Optimal beta:     " << setprecision(6) << beta << "\n"
        << pfx << "Optimal wavelength: " << setprecision(3) << wavelength << "\n";

    for(int i=0; i<nlayers; i++) {
        cout << pfx << " Layer " << setw(2) << i+1 << ": n="
            << setprecision(4) << vn[i] << setprecision(2) <<" t=" << vt[i];
        cout << "\n";
    }
    cout << flush;
    cout.unsetf(ios::fixed|ios::left);
    cout.width(0);
}

double MLADescription::CalculateTotalThickness()
{
    totalThickness=0;
    for(int i=0; i<nlayers; i++)
        totalThickness+=vt[i];
    return totalThickness;
}

double MLADescription::GetMinimumThickness() const
{
    if( nlayers==0 ) return 0.;

    double t=vt[0];

    for(size_t i=0; i<nlayers; i++)
        if( t>vt[i] ) t=vt[i];

    return t;
}

#include <gsl/gsl_sf_expint.h> //for exponential integral

MLAResult& MLADescription::Calculate(Spectrum& eff,double ps,double b)
{
    static const int Nsp=20;
    static const double eulergamma=0.5772156649;
    //normalizing coefficient of Cerenkov emission intensity assuming wavelength in nm and path in mm
    static const double K=2*M_PI/137.036/1e-6;
    static double rn2[Nsp][100];
    assert(nlayers<100);

    result.valid=false;

    double wl1, wl2;
    eff.GetRange(wl1,wl2);

    if( b==0.0 ) b=beta;

    double b2=b*b;

    //Find minimum and maximum radii of the ring
    double rnmin=b*AerogelRefIndex(vn[0],400.,wl2), rnmax=b*AerogelRefIndex(vn[0],400.,wl1);
    rn2[0][0]=rnmin*rnmin;
    rn2[1][0]=rnmax*rnmax;
    double Rmin=0, Rmax=0;
    if( rnmin>1.0 )
        Rmin=t0*sqrt((rn2[0][0]-1)/(b2-rn2[0][0]+1));
    if( rnmax>1.0 )
        Rmax=vt[0]*sqrt(rn2[1][0]-1)+t0*sqrt((rn2[1][0]-1)/(b2-rn2[1][0]+1));

    for(int l=1; l<nlayers; l++) {
        rnmin=b*AerogelRefIndex(vn[l],400.,wl2), rnmax=b*AerogelRefIndex(vn[l],400.,wl1);
        rn2[0][l]=rnmin*rnmin;
        rn2[1][l]=rnmax*rnmax;
        if( rnmax<1.0 ) continue;
        double R1=0, R2=t0*sqrt((rn2[1][l]-1)/(b2-rn2[1][l]+1));
        if( rnmin>1.0 ) R1=t0*sqrt((rn2[0][l]-1)/(b2-rn2[0][l]+1));
        for(int pl=0; pl<l; pl++) {
            if( rnmin>1.0 )
                R1+=vt[pl]*sqrt((rn2[0][l]-1)/(rn2[0][pl]-rn2[0][l]+1));
            R2+=vt[pl]*sqrt((rn2[1][l]-1)/(rn2[1][pl]-rn2[1][l]+1));
        }
        R2+=vt[l]*sqrt(rn2[1][l]-1);
        if( R1<Rmin ) Rmin=R1;
        if( R2>Rmax ) Rmax=R2;
    }

    if( Rmax==0 ) return result;

    if( !finite(Rmin) || !finite(Rmax) )
        cout<<"Rmin="<<Rmin<<" Rmax="<<Rmax<<endl;

    double Rstep=(Rmax-Rmin)/(Nr-1);

    //Preparing calculation to save CPU time
    double rn;
    double wlstep=(wl2-wl1)/(Nsp-1);
    double wl[Nsp];
    double lsc[Nsp];
    for(int iwl=0; iwl<Nsp; iwl++) {
        wl[iwl]=wl1+iwl*wlstep;
        lsc[iwl]=scatteringLength*pow(wl[iwl]/400.,4);
        for(int l=0; l<nlayers; l++) {
            rn=b*AerogelRefIndex(vn[l],400.,wl[iwl]);
            rn2[iwl][l]=rn*rn;
        }
    }

    double Npe=0;
    double Rmean=0, R2mean=0;
    double R, S, Swl, path, dR, tanc, x, att;

    for(int i=0; i<Nr; i++) {
        double R=Rmin+i*Rstep;
        double S=0;
        double Swl;
        for(int iwl=0; iwl<Nsp; iwl++) {
            Swl=0;
            for(int l=0; l<nlayers; l++) {
                if( rn2[iwl][l]<1.0 ) continue; //underthreshold velocity
                dR=t0*sqrt((rn2[iwl][l]-1)/(b2-rn2[iwl][l]+1)); //shift of the ring
                path=0; //path of light in the forward layers
                for(int pl=0; pl<l; pl++) {
                    dR+=vt[pl]*sqrt((rn2[iwl][l]-1)/(rn2[iwl][pl]-rn2[iwl][l]+1));
                    path+=vt[pl]/sqrt(1-(rn2[iwl][l]-1)/rn2[iwl][pl]);
                }

                tanc=sqrt(rn2[iwl][l]-1);

                x=(R-dR)/tanc; //position of photon emission in the current layer

                if( x<0 || x>vt[l] ) //no contribution into the current point from this layer
                    continue;

                path+=x*sqrt(rn2[iwl][l]);

                if( scatteringLength>0 )
                    att=exp(-path/lsc[iwl]);
                else
                    att=1.0;

                Swl+=(1-1/rn2[iwl][l])/tanc*att;

            } //loop on layers

            S+=K*0.01*eff.Evaluate(wl[iwl])*Swl*wlstep/wl[iwl]/wl[iwl];

        } //loop on wavelength

        result.r[i]=R;
        result.s[i]=S;

        if( S!=0 ) {
            Npe+=S*Rstep;
            Rmean+=R*S*Rstep;
            R2mean+=R*R*S*Rstep;
        }

    } //loop on radius

    if( !finite(Npe) || Npe<=0 ) {
        cout<<"Npe="<<Npe<<endl;
        return result;
    }

    Rmean/=Npe;
    R2mean/=Npe;

    result.npe=Npe;
    result.radius=Rmean;
    result.sigma1=sqrt(R2mean-Rmean*Rmean);
    if( ps!=0 ) result.sigma1=sqrt(result.sigma1*result.sigma1+ps*ps/12);
    //Exact calculation of the error per track
    double factor=(gsl_sf_expint_Ei(Npe)-eulergamma-log(Npe))/(exp(Npe)-1);
    result.sigma_t=result.sigma1*sqrt(factor);
    result.rmin=Rmin;
    result.rmax=Rmax;
    result.valid=true;

    return result;
}

