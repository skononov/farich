#include "MLADescription.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

#include "gsl/gsl_poly.h"      //for polynomial functions
#include "gsl/gsl_sf_expint.h" //for exponential integral

#include "Spectrum.h"

#include "TDecompSVD.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::setprecision;
using std::setw;

double MLADescription::tolerance = 1e-6;

ROOT::Math::Minimizer *MLADescription::minimizer = nullptr;

void MLADescription::InitializeMinimizer()
{
    if (!minimizer) {
        minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetMaxFunctionCalls(100000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.01);
        minimizer->SetPrintLevel(1);
    } else {
        minimizer->Clear();
    }
}

void MLADescription::Resolution::resize(int nl)
{
    r.resize(Nr);
    s.resize(Nr);
    data.resize(Nr * Nwl * nl);
}

MLADescription::MLADescription(double d, double b, double wl)
    : nlayers(0), t0(d), beta(b), wavelength(wl), scatteringLength(0), pixelSize(0.), optimization(0), chromaticity(true), result(),
      nPolCoef(0)
{
}

MLADescription::MLADescription(const MLADescription &mla)
    : nlayers(mla.nlayers), t0(mla.t0), beta(mla.beta), wavelength(mla.wavelength),
      scatteringLength(mla.scatteringLength), pixelSize(mla.pixelSize), optimization(mla.optimization), 
      chromaticity(mla.chromaticity), vn(mla.vn), vt(mla.vt), result(mla.result), nPolCoef(mla.nPolCoef), polCoef(mla.polCoef)
{
}

MLADescription::~MLADescription() {}

void MLADescription::Clear()
{
    rn.clear();
    rt.clear();
    vn.clear();
    vt.clear();
    polCoef.clear();
    nPolCoef = 0;
    nlayers = 0;
    optimization = 0;
    result.valid = false;
    result.r.clear();
    result.s.clear();
    result.data.clear();
}

void MLADescription::GoToAbs()
{
    nlayers = rn.size() - 1;
    vn.resize(nlayers);
    vt.resize(nlayers);
    for (int i = 0; i < nlayers; i++) {
        vn[i] = AerogelRefIndex(rn[i + 1] / beta, wavelength, 400.);
        vt[i] = rt[i + 1] * t0;
    }
}

int MLADescription::MakeAlayer()
{
    int num = rn.size();

    rn.push_back(rn.back()); // zero approximation

    int nit = 0;
    double S, disc;
    do { // iterate
        S = 0;
        for (int i = 0; i < num; i++) {
            S += rt[i] * Skl(num, i);
        }
        disc = T10 - S * Tkl(num, num);
        S = T10 / S;
        rn[num] = sqrt(1 + S * S);
        nit++;
        if (nit > 500) {
            rn.pop_back(); // reset this calculation
            return 1;      // too many loops
        }
    } while (fabs(disc) > T10 * tolerance);

    rt.push_back(rt[1] * Tkl(1, 1) / Tkl(num, num)); // derive thickness of the layer
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

double MLADescription::AerogelRefIndex(double n, double wl1, double wl2) const
{
    // Parameters of the one-pole Sellmeier formula fit to Novosibirsk aerogel n=1.03 data.
    // T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
    // Eur. Phys. J. C 52 (2007) 759-764
    static const double LHCb_a0 = 0.05639, LHCb_wl0sqr = 83.22 * 83.22;

    if (chromaticity) {
        double LHCb_RI2m1ref = LHCb_a0 / (1 - LHCb_wl0sqr / (wl1 * wl1)); //(n**2-1) of LHCb aerogel at wl1
        double ri2m1_lhcb = LHCb_a0 / (1 - LHCb_wl0sqr / (wl2 * wl2));    //(n**2-1) of LHCb aerogel at wl2
        return sqrt(1 + (n * n - 1) / LHCb_RI2m1ref * ri2m1_lhcb);
    } else {
        return n;
    }
}

void MLADescription::AddAlayer(double ri, double t)
{
    optimization = 0;
    nlayers++;
    vn.push_back(ri);
    vt.push_back(t);
}

bool MLADescription::MakeLayers(int N, double n1, double t1)
{ // fixed number of layers. Input parameters: beta, t0, n1 (at 400nm), t1, N.
    if (N < 1 || N > Nlmax) {
        cerr << "MLADescription::MakeLayers(): Error: invalid number of layers " << N << ". Should be between 1 and "
             << Nlmax << "." << endl;
        return false;
    }
    if (n1 < 1.) {
        cerr << "MLADescription::MakeLayers(): Error: refractive index " << n1 << " is less than 1.0" << endl;
        return false;
    }
    if (t1 < 0.) {
        cerr << "MLADescription::MakeLayers(): Error: i thickness value " << t1 << endl;
        return false;
    }

    Clear();

    rn.reserve(N + 1);
    rt.reserve(N + 1);

    if (wavelength != 400.)
        n1 = AerogelRefIndex(n1, 400., wavelength);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1 * beta);
    rt.push_back(t1 / t0);

    T10 = Tkl(1, 0);

    int i = 1;

    while (i < N) { // loop on layer number
        // derive index of the i-th layer
        if (MakeAlayer() != 0)
            return false;
        i++;
    }
    
    GoToAbs();

    return true;
}

bool MLADescription::MakeGabarit(double G, double n1, double rt1)
{ // fixed total thickness of radiator. Input parameters: beta, t0, G, n1 (at 400nm), rt1.
    if (n1 < 1.) {
        cerr << "MLADescription::MakeGabarit(): Error: refractive index " << n1 << " is less than 1.0" << endl;
        return false;
    }
    if (rt1 < 0.) {
        cerr << "MLADescription::MakeGabarit(): Error: invalid relative thickness value " << rt1 << endl;
        return false;
    }
    if (G < t0) {
        cerr << "MLADescription::MakeGabarit(): Error: dimension " << G << " is less than proximity distance " << t0
             << endl;
        return false;
    }

    double t1 = rt1 * t0; // thickness of the first layer

    double Tmax = G - t0; // maximum thickness of radiator (may be unachievable)

    Clear();

    if (wavelength != 400.)
        n1 = AerogelRefIndex(n1, 400., wavelength);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1 * beta);
    rt.push_back(rt1);

    T10 = Tkl(1, 0);

    double T = 0, delta = t1;
    int i = 0;

    do { // loop until required thickness gained
        T += delta;
        if (MakeAlayer() != 0)
            return false;           // something wrong with calculation
        delta = t0 * rt.back(); // thickness of a latest layer evaluated
        i++;
    } while (T + delta < Tmax);

    rt.pop_back(); // remove the last layer
    rn.pop_back(); // remove the last layer
    i--;

    GoToAbs();

    return true;
}

bool MLADescription::MakeFixed(int N, double G, double n1)
{ // fixed number of layers and total thickness of radiator.
    // Input parameters: beta, t0, N, G, n1 (at 400nm)
    if (N < 1 || N > Nlmax) {
        cerr << "MLADescription::MakeFixed(): Error: invalid number of layers " << N << ". Should be between 1 and "
             << Nlmax << "." << endl;
        return false;
    }
    if (n1 < 1.) {
        cerr << "MLADescription::MakeFixed(): Error: refractive index " << n1 << " is less than 1.0" << endl;
        return false;
    }
    if (G < t0) {
        cerr << "MLADescription::MakeFixed(): Error: dimension " << G << " is less than proximity distance " << t0
             << endl;
        return false;
    }

    Clear();

    if (N <= 0)
        return 0;

    if (N == 1) { // single layer case
        nlayers = 1;
        double t1 = G - t0;
        vn.push_back(n1);
        vt.push_back(t1);
        return true;
    }

    if (wavelength != 400.)
        n1 = AerogelRefIndex(n1, 400., wavelength);

    rn.reserve(N + 1);
    rt.reserve(N + 1);

    rn.push_back(beta);
    rt.push_back(1.0);
    rn.push_back(n1 * beta);
    T10 = Tkl(1, 0);

    double Tgoal = G - t0; // total thickness of radiator to attain
    double T = 0;

    double t1 = Tgoal / N; // estimation of the 1-st layer thickness

    int nit = 0;

    while (fabs(T - Tgoal) > tolerance * Tgoal) { // loop until required accuracy of thickness gained
        rn.resize(2);
        rt.resize(1);

        rt.push_back(t1 / t0);
        T = t1;
        int i = 1;
        for (; i < N; i++) {
            // derive the i-th layer
            if (MakeAlayer() != 0)
                break;
            T += t0 * rt[i + 1];
        }
        if (i < N)
            return false; // layers calculation failed

        t1 += (Tgoal - T) * t1 / T;

        nit++;
        if (nit > 500)
            return false; // iterations diverging
    }

    GoToAbs();

    return true;
}

void MLADescription::Print(const char *pfx) const
{
    if (!nlayers) {
        cout << "Radiator not defined yet!" << endl;
        return;
    }

    cout.setf(std::ios::fixed | std::ios::left);
    cout << pfx << "Number of layers: " << std::setprecision(3) << nlayers << "\n"
         << pfx << "Proximity distance:  " << t0 << "\n"
         << pfx << "Total thickness:  " << GetTotalThickness() << "\n"
         << pfx << "Optimal beta:     " << setprecision(6) << beta << "\n"
         << pfx << "Optimal wavelength: " << setprecision(3) << wavelength << "\n";

    if (optimization == POL_OPTIMIZATION) {
        cout.unsetf(std::ios::fixed);
        cout.setf(std::ios::scientific);
        cout << pfx << "Polynomial description: n(x1)";
        for (int i = 1; i <= nPolCoef; i++) {
            cout << " + C" << i << "*(x-x1)";
            if (i > 1)
                cout << "^" << i;
        }
        cout << "\n";
        for (int i = 1; i <= nPolCoef; i++)
            cout << pfx << " C" << i << "=" << setprecision(4) << polCoef[i-1] << "\n";
        cout.setf(std::ios::fixed);
        cout.unsetf(std::ios::scientific);
    }
    cout << pfx << "Layer thicknesses and refractive indices (20 layers at most):\n";
    for (int l = 0; l < std::min(20, nlayers); l++) {
        cout << pfx << " Layer" << setw(3) << l + 1 << ": n=" << setprecision(4) << vn[l] << setprecision(2)
             << " t=" << vt[l] << "\n";
    }
    if (nlayers > 20)
        cout << pfx << " ......\n";
    cout << flush;
    cout.unsetf(std::ios::fixed | std::ios::left);
    cout.precision(6);
    cout.width(0);
}

double MLADescription::GetMaxSensitivityWL(double T) const
{
    double wl1, wl2;
    double wlmax = 400, smax = 0;
    double wl, rn, Latt, attFactor, s;

    pdEff.GetRange(wl1, wl2);
    double wlstep = (wl2 - wl1) / Nwl;

    for (int i = 0; i < Nwl; i++) {
        wl = wl1 + (i + 0.5) * wlstep;
        rn = AerogelRefIndex(1.05, 400., wl) * beta;
        if (rn <= 1.0)
            continue;
        attFactor = 1.;
        if (T > 0.) {
            Latt = scatteringLength * pow(wl / 400., 4);
            if (!absLength.IsEmpty()) Latt = 1 / (1/absLength(wl) + 1/Latt);
            attFactor = Latt * (1 - exp(- T * rn / Latt));
        }
        s = (1 - 1 / (rn * rn)) * attFactor * pdEff.Evaluate(wl) / rn;
        if (smax < s) {
            smax = s;
            wlmax = wl;
        }
    }

    return wlmax;
}

MLADescription::Resolution &MLADescription::Calculate(double b, bool storeData, bool fresnelOn)
{
    static const double eulergamma = 0.5772156649;
    // normalizing coefficient of Cerenkov emission intensity assuming wavelength in nm and path in mm
    static const double K = 2 * M_PI / 137.036 / 1e-6;
    static double rn2[Nwl][Nlmax], tanp[Nwl][Nlmax][Nlmax], cosp[Nwl][Nlmax][Nlmax], tana[Nwl][Nlmax];
    static double wl[Nwl];
    static double lsc[Nwl];

    result.valid = false;

    double wl1, wl2;
    pdEff.GetRange(wl1, wl2);

    if (b <= 0.)
        b = beta;
    result.beta = b;

    double b2 = b * b;

    // Find minimum and maximum radii of the ring
    double rnmin = b * AerogelRefIndex(vn[0], 400., wl2), rnmax = b * AerogelRefIndex(vn[0], 400., wl1);
    rn2[0][0] = rnmin * rnmin;
    rn2[1][0] = rnmax * rnmax;
    double Rmin = 0, Rmax = 0;
    if (rnmin > 1.0)
        Rmin = t0 * sqrt((rn2[0][0] - 1) / (b2 - rn2[0][0] + 1));
    if (rnmax > 1.0)
        Rmax = vt[0] * sqrt(rn2[1][0] - 1) + t0 * sqrt((rn2[1][0] - 1) / (b2 - rn2[1][0] + 1));

    for (int l = 1; l < nlayers; l++) {
        rnmin = b * AerogelRefIndex(vn[l], 400., wl2), rnmax = b * AerogelRefIndex(vn[l], 400., wl1);
        rn2[0][l] = rnmin * rnmin;
        rn2[1][l] = rnmax * rnmax;
        if (rnmax < 1.0)
            continue;
        double R1 = 0, R2 = t0 * sqrt((rn2[1][l] - 1) / (b2 - rn2[1][l] + 1));
        if (rnmin > 1.0)
            R1 = t0 * sqrt((rn2[0][l] - 1) / (b2 - rn2[0][l] + 1));
        for (int pl = 0; pl < l; pl++) {
            if (rnmin > 1.0)
                R1 += vt[pl] * sqrt((rn2[0][l] - 1) / (rn2[0][pl] - rn2[0][l] + 1));
            R2 += vt[pl] * sqrt((rn2[1][l] - 1) / (rn2[1][pl] - rn2[1][l] + 1));
        }
        R2 += vt[l] * sqrt(rn2[1][l] - 1);
        if (R1 < Rmin)
            Rmin = R1;
        if (R2 > Rmax)
            Rmax = R2;
    }

    if (Rmax == 0)
        return result;

    assert(finite(Rmin));
    assert(finite(Rmax));
    assert(Rmin < Rmax);

    if (storeData)
        result.resize(nlayers);

    double Rstep = (Rmax - Rmin) / Nr;

    // Preparing calculation to save CPU time
    double wlstep = (wl2 - wl1) / Nwl;

#pragma omp parallel for
    for (int iwl = 0; iwl < Nwl; iwl++) {
        wl[iwl] = wl1 + (iwl + 0.5) * wlstep;
        lsc[iwl] = scatteringLength * pow(wl[iwl] / 400., 4);
        for (int l = 0; l < nlayers; l++) {
            double rn = b * AerogelRefIndex(vn[l], 400., wl[iwl]), rnsq = rn * rn;
            rn2[iwl][l] = rnsq;
            if (rnsq <= 1.0) { // underthreshold velocity
                tanp[iwl][l][l] = 0;
                tana[iwl][l] = 0;
                continue;
            }
            tana[iwl][l] = sqrt((rnsq - 1) / (b2 - rnsq + 1));
            for (int pl = 0; pl <= l; pl++) {
                tanp[iwl][l][pl] = sqrt((rnsq - 1) / (rn2[iwl][pl] - rnsq + 1));
                cosp[iwl][l][pl] = sqrt(1 - (rnsq - 1) / rn2[iwl][pl]);
            }
        }
    }

    double Npe = 0, Rmean = 0, R2mean = 0;

#pragma omp parallel for reduction(+ : Npe, Rmean, R2mean)
    for (int i = 0; i < Nr; i++) { // loop on radius
        double R = Rmin + i * Rstep, Rc = R + 0.5 * Rstep;
        double S = 0;
        for (int iwl = 0; iwl < Nwl; iwl++) { // loop on wavelength
            double Swl = 0, X0 = 0;
            for (int l = 0; l < nlayers; X0 += vt[l], l++) { // loop on layers
                if (storeData)
                    result.data[l * Nwl * Nr + iwl * Nr + i] = {l, (float)Rc, (float)wl[iwl], 0., 0.};
                if (rn2[iwl][l] <= 1.0) // underthreshold velocity
                    continue;
                double dR = t0 * tana[iwl][l]; // shift of the ring in air
                double path = 0; // path of light in the forward layers
                double att = 1.0;
                if (fresnelOn) 
                    att *= 1. - pow((tanp[iwl][l][0]-tana[iwl][l]) / (tanp[iwl][l][0]+tana[iwl][l]) * 
                                    (1-tanp[iwl][l][0]*tana[iwl][l]) / (1+tanp[iwl][l][0]*tana[iwl][l]),2);
                for (int pl = 0; pl < l; pl++) {
                    dR += vt[pl] * tanp[iwl][l][pl];
                    path += vt[pl] / cosp[iwl][l][pl];
                    if (fresnelOn)
                        att *= 1. - pow((tanp[iwl][l][pl]-tanp[iwl][l][pl+1]) / (tanp[iwl][l][pl]+tanp[iwl][l][pl+1]) *
                                        (1-tanp[iwl][l][pl]*tanp[iwl][l][pl+1]) / (1+tanp[iwl][l][pl]*tanp[iwl][l][pl+1]),2);
                }

                double x0 = (R - dR) / tanp[iwl][l][l];   // position of photon emission for the first radius point
                double x1 = x0 + Rstep / tanp[iwl][l][l]; // position of photon emission for the second radius point

                if (x0 < 0. && x1 < 0. ||
                    x0 > vt[l] && x1 > vt[l]) // no contribution into the current point from this layer
                    continue;

                x0 = x0 < 0. ? 0. : x0;
                x1 = x1 > vt[l] ? vt[l] : x1;

                path += 0.5 * (x0 + x1) / cosp[iwl][l][l];

                if (scatteringLength > 0)
                    att *= exp(-path / lsc[iwl]);
                if (!absLength.IsEmpty())
                    att *= exp(-path / absLength(wl[iwl]));

                double dSwl = (1 - 1 / rn2[iwl][l]) * att * (x1 - x0);
                if (storeData)
                    result.data[l * Nwl * Nr + iwl * Nr + i] = {
                        l, (float)Rc, (float)wl[iwl], (float)(X0 + 0.5 * (x0 + x1)),
                        (float)(K * 0.01 * pdEff(wl[iwl]) * dSwl * wlstep / wl[iwl] / wl[iwl])};

                Swl += dSwl;
            } // loop on layers

            S += K * 0.01 * pdEff.Evaluate(wl[iwl]) * Swl * wlstep / wl[iwl] / wl[iwl];
        } // loop on wavelength

        if (storeData) {
            result.r[i] = Rc;
            result.s[i] = S / Rstep;
        }

        Npe += S;
        Rmean += Rc * S;
        R2mean += Rc * Rc * S;
    } // loop on radius

    if (!finite(Npe) || Npe <= 0) {
        cout << "MLADescription::Calculate(): Npe=" << Npe << endl;
        Print("!!!");
        return result;
    }

    Rmean /= Npe;
    R2mean /= Npe;

    result.npe = Npe;
    result.rmin = Rmin;
    result.rmax = Rmax;
    result.rstep = Rstep;

    result.radius = Rmean;

    result.sigma1 = sqrt(R2mean - Rmean * Rmean);
    result.sigma1_px = sqrt(R2mean - Rmean * Rmean + pixelSize * pixelSize / 12);

    double D = t0 + 0.5 * GetTotalThickness(), fmm_mrad = 1e3 * D / (D * D + Rmean * Rmean);

    result.sigma1_ang = result.sigma1 * fmm_mrad;
    result.sigma1_ang_px = result.sigma1_px * fmm_mrad;

    // Exact calculation of the error per track
    double factor = sqrt((gsl_sf_expint_Ei(Npe) - eulergamma - log(Npe)) / (exp(Npe) - 1));

    result.sigma_t = result.sigma1 * factor;
    result.sigma_t_px = result.sigma1_px * factor;

    result.sigma_t_ang = result.sigma_t * fmm_mrad;
    result.sigma_t_ang_px = result.sigma_t_px * fmm_mrad;

    result.valid = true;

    return result;
}

void MLADescription::ApplyNTParameterization(const double *p)
{
    TVectorD pt(nlayers);
    for (int l = 0; l < nlayers; l++)
        pt[l] = p[2 * l + 1];

    TVectorD t = V * pt;

    for (int l = 0; l < nlayers; l++) {
        vn[l] = p[2 * l];
        vt[l] = std::max(0., t[l]);
    }
}

double MLADescription::EvalResolutionNT(const double *p)
{
    ApplyNTParameterization(p);

    Calculate();

    if (!result.valid)
        return NAN;

    return result.sigma_t_ang;
}

void MLADescription::ApplyPolParameterization(const double *p)
{
    double dx = 0.;
    std::copy(p, p+nPolCoef, polCoef.begin());
    if (nPolCoef == 1) {
        for (int l = 1; l < nlayers; l++) {
            dx += 0.5 * (vt[l - 1] + vt[l]);
            vn[l] = std::max(1., vn[0] + dx * p[0]);
        }
    } else {
        for (int l = 1; l < nlayers; l++) {
            dx += 0.5 * (vt[l - 1] + vt[l]);
            // polynomial with fixed value of refractive index in the middle of the first layer
            vn[l] = std::max(1., vn[0] + dx * gsl_poly_eval(p, nPolCoef, dx)); 
        }
    }
}

double MLADescription::EvalResolutionPol(const double *p)
{
    ApplyPolParameterization(p);

    Calculate();

    if (!result.valid)
        return NAN;

    return result.sigma_t_ang;
}

bool MLADescription::OptimizeNT(int N, double G, double n1)
{
    // Make focusing radiator via analytical calculations as the first approximaiton
    if (!MakeFixed(N, G, n1)) return false;
    
    // Skip optimization for a single layer
    if (N==1) return true;

    cout << "Optimize radiator with NT parameterization" << endl;

    // Prepare calculations evaluating V and Vtr transormation matrices
    TVectorD t(nlayers), pt(nlayers);
    for (int i = 0; i < nlayers; i++)
        t[i] = vt[i];

    TMatrixD A(nlayers, nlayers);
    TVectorD ones(nlayers);
    ones = 1;
    TMatrixDRow(A, 0) = ones;

    TDecompSVD D(A);
    D.Decompose();

    V.ResizeTo(A);
    Vtr.ResizeTo(A);

    V = D.GetV();

    Vtr = V;
    Vtr.T();

    pt = Vtr * t;

    InitializeMinimizer();

    optimization = NT_OPTIMIZATION;

    // Define errors for 10% resolution variation based on current radiator resolution
    Calculate();
    minimizer->SetErrorDef(0.1*result.sigma_t_ang);

    ROOT::Math::Functor fcn(this, &MLADescription::EvalResolutionNT, 2 * nlayers);

    minimizer->SetFunction(fcn);
    // Set the free variables
    char name[20];
    for (int l = 0; l < nlayers; l++) {
        sprintf(name, "n%d", l + 1);
        minimizer->SetLimitedVariable(2 * l, name, vn[l], 0.005, 1.0, 1.2);
        sprintf(name, "pt%d", l + 1);
        minimizer->SetVariable(2 * l + 1, name, pt[l], 0.5);
    }
    minimizer->FixVariable(0);
    minimizer->FixVariable(1);

    minimizer->Minimize();

    if (minimizer->Status() > 1) {
        cout << ">> Failed to optimize. Status = " << minimizer->Status() << endl;
        return false;
    }

    ApplyNTParameterization(minimizer->X());

    return true;
}

bool MLADescription::OptimizePol(int N, int npol, double G, double nmax, bool sameThick)
{
    if (N < 1 || N > Nlmax) {
        cerr << "MLADescription::OptimizePol(): Error: invalid number of layers " << N << ". Should be between 1 and "
             << Nlmax << "." << endl;
        return false;
    }
    if (npol < 1) {
        cerr << "MLADescription::OptimizePol(): Error: polynomial degree " << npol
             << " is too little. Should be at least 1." << endl;
        return false;
    }
    if (G < t0) {
        cerr << "MLADescription::OptimizePol(): Error: dimension " << G << " is less than proximity distance " << t0
             << endl;
        return false;
    }
    if (nmax < 1.) {
        cerr << "MLADescription::OptimizePol(): Error: refractive index " << nmax << " is less than 1.0" << endl;
        return false;
    }

    // Estimate initial values of parameters by linear fitting of the fast optimized radiator
    if (!MakeFixed(N, G, nmax)) return false;
    
    // Skip optimization for a single layer
    if (N == 1) return true;    

    TVectorD C(npol);
    if (N > 2) {
        TVectorD b(N-1);
        TMatrixD A(N-1,npol);
        float dx = 0.;
        for(int l=1; l<N; l++) {
            dx += 0.5 * (vt[l-1] + vt[l]);
            b[l-1] = (vn[l]-nmax)/dx;
            for(int k=0; k<npol; k++)
                A[l-1][k] = pow(dx,k);
        }
        
        TDecompSVD svd(A);
        Bool_t ok=kTRUE;
    
        C = svd.Solve(b,ok);
    
        if (!ok) {
            cerr << "MLADescription::OptimizePol(): Error: Can not determine initial parameter values" << endl;
            return false;
        }
    } else { // N==2
        // Gradient estimation (negative)
        C[0] = 2. * (vn[1]-vn[0]) / (vt[0]+vt[1]);
    }
    
    cout << "Optimize radiator with polynomial parameterization" << endl;

    double T = G - t0; // total radiator thickness

    // Make layers of the same thickness if sameThick is true
    if (sameThick) {
        Clear();
        nlayers = N;
        vn.resize(nlayers, nmax);
        vt.resize(nlayers, T / nlayers);
    }

    nPolCoef = npol;
    polCoef.resize(nPolCoef);

    InitializeMinimizer();

    optimization = POL_OPTIMIZATION;

    ApplyPolParameterization(C.GetMatrixArray());

    // Define errors for 10% resolution variation based on current radiator resolution
    Calculate();
    minimizer->SetErrorDef(0.1*result.sigma_t_ang);

    ROOT::Math::Functor fcn(this, &MLADescription::EvalResolutionPol, nPolCoef);

    minimizer->SetFunction(fcn);
    // Set the free variables
    char name[20];
    for (int i = 0; i < nPolCoef; i++) {
        sprintf(name, "C%d", i+1);
        double initStep = 0.01*std::max(fabs(polCoef[i]), fabs(polCoef[0])/pow(T,i)); // should be small but not zero
        minimizer->SetVariable(i, name, polCoef[i], initStep);
    }

    minimizer->Minimize();

    if (minimizer->Status() > 1) {
        cout << ">> Failed to optimize. Status = " << minimizer->Status() << endl;
        return false;
    }

    ApplyPolParameterization(minimizer->X());

    return true;
}
