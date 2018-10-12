#include "MLADescription.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

#include "gsl/gsl_poly.h" //for polynomial functions
#include "gsl/gsl_sf_expint.h" //for exponential integral

#include "Spectrum.h"

#include "TDecompSVD.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#ifdef _OPENMP
# include <omp.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setprecision;
using std::setw;

double MLADescription::tolerance = 1e-6;

ROOT::Math::Minimizer* MLADescription::minimizer = nullptr;

void MLADescription::InitializeMinimizer()
{
    if( !minimizer ) {
        minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetMaxFunctionCalls(1000000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.001);
        minimizer->SetPrintLevel(1);
    } else {
        minimizer->Clear();
    }
}

void MLADescription::Resolution::reserve(int nl)
{
    r.clear();
    s.clear();
    data.clear();

    r.reserve(Nr);
    s.reserve(Nr);
    data.reserve(Nr*Nwl*nl);
}

MLADescription::MLADescription(double d, double b, double wl) :
    nlayers(0),
    t0(d),
    beta(b),
    wavelength(wl),
    scatteringLength(0),
    pixelSize(0.),
    optimization(0),
    result(),
    nPolCoef(0)
{
}

MLADescription::MLADescription(const MLADescription& mla) :
    nlayers(mla.nlayers),
    t0(mla.t0),
    beta(mla.beta),
    wavelength(mla.wavelength),
    scatteringLength(mla.scatteringLength),
    pixelSize(mla.pixelSize),
    optimization(mla.optimization),
    vn(mla.vn),
    vt(mla.vt),
    result(mla.result),
    nPolCoef(mla.nPolCoef),
    polCoef(mla.polCoef)
{
}

MLADescription::~MLADescription()
{}

void MLADescription::Clear()
{
    rn.clear();
    rt.clear();
    vn.clear();
    vt.clear();
    polCoef.clear();
    nPolCoef=0;
    nlayers=0;
    optimization=0;
    result.valid=false;
    result.r.clear();
    result.s.clear();
    result.data.clear();
}

void MLADescription::GoToAbs()
{
    nlayers=rn.size()-1;
    vn.resize(nlayers);
    vt.resize(nlayers);
    for(int i=0; i<nlayers; i++) {
        vn[i]=AerogelRefIndex(rn[i+1]/beta,wavelength,400.);
        vt[i]=rt[i+1]*t0;
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
    optimization=0;
    nlayers++;
    vn.push_back(ri);
    vt.push_back(t);
}

bool MLADescription::MakeLayers(int N, double n1, double t1)
{ //fixed number of layers. Input parameters: beta, t0, n1 (at 400nm), t1, N.
    if( N<1 || N>Nlmax ) {
        cerr<<"MLADescription::MakeLayers(): Error: invalid number of layers "<<N<<". Should be between 1 and "<<Nlmax<<"."<<endl;
        return false;
    }
    if( n1<1. ) {
        cerr<<"MLADescription::MakeLayers(): Error: refractive index "<<n1<<" is less than 1.0"<<endl;
        return false;
    }
    if( t1<0. ) {
        cerr<<"MLADescription::MakeLayers(): Error: i thickness value "<<t1<<endl;
        return false;
    }

    Clear();

    rn.reserve(N+1);
    rt.reserve(N+1);

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

    return true;
}

bool MLADescription::MakeGabarit(double G, double n1, double rt1)
{ //fixed total thickness of radiator. Input parameters: beta, t0, G, n1 (at 400nm), rt1.
    if( n1<1. ) {
        cerr<<"MLADescription::MakeGabarit(): Error: refractive index "<<n1<<" is less than 1.0"<<endl;
        return false;
    }
    if( rt1<0. ) {
        cerr<<"MLADescription::MakeGabarit(): Error: invalid relative thickness value "<<rt1<<endl;
        return false;
    }
    if( G<t0 ) {
        cerr<<"MLADescription::MakeGabarit(): Error: dimension "<<G<<" is less than proximity distance "<<t0<<endl;
        return false;
    }
    
    double t1 = rt1*t0; //thickness of the first layer

    double Tmax=G-t0; //maximum thickness of radiator (may be unachievable)

    Clear();

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

    return true;
}

bool MLADescription::MakeFixed(int N, double G, double n1)
{ //fixed number of layers and total thickness of radiator.
    //Input parameters: beta, t0, N, G, n1 (at 400nm)
    if( N<1 || N>Nlmax ) {
        cerr<<"MLADescription::MakeFixed(): Error: invalid number of layers "<<N<<". Should be between 1 and "<<Nlmax<<"."<<endl;
        return false;
    }
    if( n1<1. ) {
        cerr<<"MLADescription::MakeFixed(): Error: refractive index "<<n1<<" is less than 1.0"<<endl;
        return false;
    }
    if( G<t0 ) {
        cerr<<"MLADescription::MakeFixed(): Error: dimension "<<G<<" is less than proximity distance "<<t0<<endl;
        return false;
    }

    Clear();

    if( N<=0 ) return 0;

    if( N==1 ) { //single layer case
        nlayers=1;
        double t1=G-t0;
        vn.push_back(n1);
        vt.push_back(t1);
        return true;
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
        if (nit>500) return false; //iterations diverging
    }

    GoToAbs();

    return true;
}

void MLADescription::Print(const char* pfx) const
{
    if( !nlayers ) {
        cout<<"Radiator not defined yet!"<<endl;
        return;
    }

    cout.setf(std::ios::fixed|std::ios::left);
    cout << pfx << "Number of layers: " << std::setprecision(3) << nlayers << "\n"
        << pfx << "Proximity distance:  " << t0 << "\n"
        << pfx << "Total thickness:  " << GetTotalThickness() << "\n"
        << pfx << "Optimal beta:     " << setprecision(6) << beta << "\n"
        << pfx << "Optimal wavelength: " << setprecision(3) << wavelength << "\n";

    if( optimization==POL_OPTIMIZATION ) {
        cout.unsetf(std::ios::fixed);
        cout.setf(std::ios::scientific);        
        cout << pfx << "Polynomial description: n(x1) + (x-x1)*(C0";
        for(int i=1; i<nPolCoef; i++) {
            cout << "+C" << i << "*x";
            if( i>1 ) cout << "^" << i;
        }            
        cout << ")\n";
        for(int i=0; i<nPolCoef; i++)
            cout << pfx << " C" << i << "=" << setprecision(4) << polCoef[i] << "\n";
        cout.setf(std::ios::fixed);
        cout.unsetf(std::ios::scientific);        
    }
    cout << pfx << "Layer thicknesses and refractive indices (20 layers at most):\n";
    for(int i=0; i<std::min(20,nlayers); i++) {
        cout << pfx << " Layer " << setw(2) << i+1 << ": n="
             << setprecision(4) << vn[i] << setprecision(2) <<" t=" << vt[i] << "\n";
    }
    if( nlayers>20 )
        cout << pfx << " ......\n";
    cout << flush;
    cout.unsetf(std::ios::fixed|std::ios::left);
    cout.precision(6);
    cout.width(0);
}

double MLADescription::GetMaxSensitivityWL(double T) const
{
    double wl1, wl2;
    double wlmax=400, smax=0;

    pdEff.GetRange(wl1,wl2);
    double wlstep=(wl2-wl1)/Nwl;    

    for(int i=0; i<Nwl; i++) {
        double wl=wl1+(i+0.5)*wlstep;
        double rn=AerogelRefIndex(vn[0],400.,wl)*beta;
        if( rn<=1.0 ) continue;
        double attFactor = 1.;
        if( T>0. ) {
            double Latt=scatteringLength*pow(wl/400.,4);
            attFactor=Latt*(1-exp(-T*rn/Latt));
        }
        double s=(1-1/(rn*rn))*attFactor*pdEff.Evaluate(wl)/rn;
        if( smax<s ) {
            smax=s;
            wlmax=wl;
        }
    }

    return wlmax;
}

MLADescription::Resolution& MLADescription::Calculate(double b, bool storeData)
{
    static const double eulergamma=0.5772156649;
    //normalizing coefficient of Cerenkov emission intensity assuming wavelength in nm and path in mm
    static const double K=2*M_PI/137.036/1e-6;
    static double rn2[Nwl][Nlmax], tanc[Nwl][Nlmax];

    result.valid=false;

    double wl1, wl2;
    pdEff.GetRange(wl1,wl2);

    if( b<=0.) b=beta;
    result.beta=b;
    
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

    if( !finite(Rmin) || !finite(Rmax) ) {
        cout<<"MLADescription::Calculate(): Rmin="<<Rmin<<" Rmax="<<Rmax<<endl;
        return result;
    }

    if( storeData )
        result.reserve(nlayers);

    double Rstep=(Rmax-Rmin)/Nr;

    //Preparing calculation to save CPU time
    double wlstep=(wl2-wl1)/Nwl;
    double wl[Nwl];
    double lsc[Nwl];
    
    #pragma omp parallel for
    for(int iwl=0; iwl<Nwl; iwl++) {
        wl[iwl]=wl1+(iwl+0.5)*wlstep;
        lsc[iwl]=scatteringLength*pow(wl[iwl]/400.,4);
        for(int l=0; l<nlayers; l++) {
            double rn=b*AerogelRefIndex(vn[l],400.,wl[iwl]);
            rn2[iwl][l]=rn*rn;
            tanc[iwl][l]=sqrt(rn*rn-1);
        }
    }

    double Npe=0, Rmean=0, R2mean=0;

    #pragma omp parallel for reduction(+:Npe,Rmean,R2mean)
    for(int i=0; i<Nr; i++) { //loop on radius
        double R=Rmin+i*Rstep, Rc=R+0.5*Rstep;
        double S=0;
        for(int iwl=0; iwl<Nwl; iwl++) { //loop on wavelength
            double Swl=0;
            for(int l=0; l<nlayers; l++) { //loop on layers
                if( rn2[iwl][l]<1.0 ) continue; //underthreshold velocity
                double dR=t0*sqrt((rn2[iwl][l]-1)/(b2-rn2[iwl][l]+1)); //shift of the ring in air
                double path=0; //path of light in the forward layers
                for(int pl=0; pl<l; pl++) {
                    dR+=vt[pl]*sqrt((rn2[iwl][l]-1)/(rn2[iwl][pl]-rn2[iwl][l]+1));
                    path+=vt[pl]/sqrt(1-(rn2[iwl][l]-1)/rn2[iwl][pl]);
                }

                double x0=(R-dR)/tanc[iwl][l]; //position of photon emission for the first radius point
                double x1=x0+Rstep/tanc[iwl][l]; //position of photon emission for the second radius point

                if( x0<0. && x1<0. || x0>vt[l] && x1>vt[l] ) //no contribution into the current point from this layer
                    continue;

                x0=x0<0.?0.:x0;
                x1=x1>vt[l]?vt[l]:x1;
                
                path+=0.5*(x0+x1)*sqrt(rn2[iwl][l]);

                double att=1.0;
                if( scatteringLength>0 )
                    att=exp(-path/lsc[iwl]);

                double dSwl = (1-1/rn2[iwl][l])*att*(x1-x0);
                if( storeData )
                    result.data.push_back({l, (float)R, (float)wl[iwl], (float)x0, (float)(K*0.01*pdEff.Evaluate(wl[iwl])*dSwl/wl[iwl]/wl[iwl])});
                
                Swl+=dSwl;
            } //loop on layers

            S+=K*0.01*pdEff.Evaluate(wl[iwl])*Swl*wlstep/wl[iwl]/wl[iwl];
        } //loop on wavelength

        if( storeData ) {
            result.r.push_back(Rc);
            result.s.push_back(S/Rstep);
        }

        Npe+=S;
        Rmean+=Rc*S;
        R2mean+=Rc*Rc*S;
    } //loop on radius

    if( !finite(Npe) || Npe<=0 ) {
        cout<<"MLADescription::Calculate(): Npe="<<Npe<<endl;
        return result;
    }

    Rmean/=Npe;
    R2mean/=Npe;

    result.npe=Npe;
    result.rmin=Rmin;
    result.rmax=Rmax;
    result.rstep=Rstep;

    result.radius=Rmean;

    result.sigma1=sqrt(R2mean-Rmean*Rmean);
    result.sigma1_px=sqrt(R2mean-Rmean*Rmean+pixelSize*pixelSize/12);

    double D=t0+0.5*GetTotalThickness(), fmm_mrad=1e3*D/(D*D+Rmean*Rmean);
        
    result.sigma1_ang=result.sigma1*fmm_mrad;
    result.sigma1_ang_px=result.sigma1_px*fmm_mrad;

    //Exact calculation of the error per track
    double factor=sqrt((gsl_sf_expint_Ei(Npe)-eulergamma-log(Npe))/(exp(Npe)-1));

    result.sigma_t=result.sigma1*factor;
    result.sigma_t_px=result.sigma1_px*factor;

    result.sigma_t_ang=result.sigma_t*fmm_mrad;
    result.sigma_t_ang_px=result.sigma_t_px*fmm_mrad;

    result.valid=true;

    return result;
}

void MLADescription::ApplyNTParameterization(const double* p)
{
    TVectorD pt(nlayers);
    for(int l=0; l<nlayers; l++)
        pt[l] = p[2*l+1];
    
    TVectorD t = V*pt;

    for(int l=0; l<nlayers; l++) {
        vn[l] = p[2*l];
        vt[l] = t[l];
    }
}

double MLADescription::EvalResolutionNT(const double* p)
{
    ApplyNTParameterization(p);
    
    Calculate();

    if( !result.valid ) return NAN;

    return result.sigma_t_ang;
}

void MLADescription::ApplyPolParameterization(const double* p)
{
    double x0 = vt[0]*0.5, x = x0;
    if( nPolCoef==1 ) {
        polCoef[0] = p[0];
        for(int l=1; l<nlayers; l++) {
            x += vt[l];
            vn[l] = std::max(1.,vn[0] + (x-x0)*p[0]);
        }
    } else {
        for(int l=1; l<nlayers; l++) {
            x += vt[l];
            vn[l] = std::max(1.,vn[0] + (x-x0)*gsl_poly_eval(p,nPolCoef,x)); //polynomial with fixed value of refractive index in the middle of the first layer
        }
    }
}

double MLADescription::EvalResolutionPol(const double* p)
{
    ApplyPolParameterization(p);

    Calculate();

    if( !result.valid ) return NAN;

    return result.sigma_t_ang;
}

bool MLADescription::OptimizeNT(int N, double G, double n1)
{
    // Make focusing radiator via analytical calculations as the first approximaiton
    MakeFixed(N, G, n1);

    cout << "Optimize radiator with NT parameterization" << endl;

    // Prepare calculations evaluating V and Vtr transormation matrices
    TVectorD t(nlayers), pt(nlayers);
    for(int i=0; i<nlayers; i++)
        t[i] = vt[i];
    
    TMatrixD A(nlayers, nlayers);
    TVectorD ones(nlayers);
    ones = 1;
    TMatrixDRow(A,0) = ones;
    
    TDecompSVD D(A);
    D.Decompose();
    
    V.ResizeTo(A);
    Vtr.ResizeTo(A);
    
    V = D.GetV();

    Vtr = V; Vtr.T();

    pt = Vtr*t;
    
    InitializeMinimizer();
    
    ROOT::Math::Functor fcn(this,&MLADescription::EvalResolutionNT,2*nlayers); 

    minimizer->SetFunction(fcn);
    // Set the free variables to be minimized
    char name[20];
    for(int l=0; l<nlayers; l++) {
        sprintf(name,"n%d",l+1);
        minimizer->SetVariable(2*l, name, vn[l], 0.005);
        sprintf(name,"pt%d",l+1);
        minimizer->SetVariable(2*l+1, name, pt[l], 0.05);
    }
    minimizer->FixVariable(0); 
    minimizer->FixVariable(1);

    minimizer->Minimize();

    if( minimizer->Status() > 1 ) {
        cout<<">> Failed to optimize. Status = " <<minimizer->Status()<<endl;
        return false;
    }

    ApplyNTParameterization(minimizer->X());

    optimization = NT_OPTIMIZATION;
    
    return true;
}

bool MLADescription::OptimizePol(int N, int npol, double G, double nmax)
{
    if( N<1 || N>Nlmax ) {
        cerr<<"MLADescription::OptimizePol(): Error: invalid number of layers "<<N<<". Should be between 1 and "<<Nlmax<<"."<<endl;
        return false;
    }
    if( npol<1 ) {
        cerr<<"MLADescription::OptimizePol(): Error: polynomial degree "<<npol<<" is too little. Should be at least 1."<<endl;
        return false;
    }
    if( G<t0 ) {
        cerr<<"MLADescription::OptimizePol(): Error: dimension "<<G<<" is less than proximity distance "<<t0<<endl;
        return false;
    }
    if( nmax<1. ) {
        cerr<<"MLADescription::OptimizePol(): Error: refractive index "<<nmax<<" is less than 1.0"<<endl;
        return false;
    }

    Clear();
    
    cout << "Optimize radiator with polynomial parameterization" << endl;

    nlayers = N;
    nPolCoef = npol;
    polCoef.resize(npol);
    double T = G-t0; //total radiator thickness
    vn.resize(nlayers, nmax);
    vt.resize(nlayers, T/nlayers);
    
    InitializeMinimizer();

    ROOT::Math::Functor fcn(this,&MLADescription::EvalResolutionPol,nPolCoef); 

    minimizer->SetFunction(fcn);
    // Set the free variables to be minimized
    char name[20];
    for(int i=0; i<nPolCoef; i++) {
        sprintf(name,"c%d",i);
        minimizer->SetVariable(i, name, 0., 1e-4/pow(T,i));
    }

    minimizer->Minimize();

    if( minimizer->Status() > 1 ) {
        cout<<">> Failed to optimize. Status = " <<minimizer->Status()<<endl;
        return false;
    }

    const double *p = minimizer->X();
    ApplyPolParameterization(p);

    //Copy final polynomial coefficients to this object
    std::copy(p,p+nPolCoef,polCoef.begin());
    
    optimization = POL_OPTIMIZATION;

    return true;
}
