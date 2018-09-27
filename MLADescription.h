#ifndef MLADescription_h
#define MLADescription_h 1

# include <vector>
# include <numeric>
# include <algorithm>
# include <string>
# include <cmath>
# include "Spectrum.h"
# include "TMatrixD.h"

namespace ROOT { namespace Math { class Minimizer; } }

/**Derived from MLRadiatorDescription class in FARICH simulation.
   Date: 27.02.2009
*/

class MLADescription {
public:
    static const int Nr=100; ///< Number of division on radius for the photon distribution calculation
    static const int Nsp=100; ///< Number of wavelength points for photon distribution calculation

    enum Optimization_t {
        NOT_OPTIMIZED = 0,
        NT_OPTIMIZATION = 1,
        POL_OPTIMIZATION = 2  
    };
    struct Resolution {
        bool valid;
        double npe;
        double radius;
        double sigma1;
        double sigma_t;
        double rmin, rmax;
        double r[Nr];
        double s[Nr];
    };

private:
    static double tolerance; ///< allowed coordinate shift calculation error

    int nlayers;              ///< number of layers in radiator
    double t0;                ///< proximity distance, mm
    double beta;              ///< particle velocity that we need to optimize for
    double wavelength;        ///< wavelength at which we need to optimize for, nm
    double scatteringLength;  ///< light scattering length at 400 nm, mm
    Spectrum pdEff;           ///< Photon detection efficiency spectrum, percent
    double pixelSize;         ///< Photon detector pixel size, mm
    int optimization;         ///< last optimization type, possible values are defined in enum
        
    //Auxilliary data structures
    std::vector<double> vn;   ///< indices of refraction at 400 nm
    std::vector<double> vt;   ///< thicknesses, mm

    std::vector<double> rn;  ///< reduced indices for calculations
    std::vector<double> rt;  ///< reduced thicknesses for calculations

    double T10;              ///< auxilary coefficient

    struct Resolution result; ///< last result of resolution calculation
        
    TMatrixD V, Vtr;         ///< V matrix from SVD decomposition of the equation: sum of layer thicknesses = totalThickness 
    int nPolCoef;            ///< Number of polynomial coefficients for refractive index profile description
    std::vector<double> polCoef; ///< Polynomic coefficients

    static ROOT::Math::Minimizer* minimizer; ///< pointer to minimizer instance

public:

    /**Set tolerance - calculation accuracy parameter*/
    static void SetTolerance(double eps) { tolerance=eps; }

    /**Constructor.
     * @param d proximity distance
     * @param b beta parameter
     */
    MLADescription(double d=1.0, double b=1.0, double wl=400);
    /**Copy constructor*/
    MLADescription(const MLADescription &mla);
    /**Destructor*/
    ~MLADescription();

    /**Print radiator configuration*/
    void Print(const char* pfx="") const;

    /**Add a layer of aerogel with refractive index @c ri, thickness @c t*/
    void AddAlayer(double ri,double t);

    /**Reset radiator*/
    void Clear();

    /**Make up focusing aerogel radiator with given number of layers.
     * @param N number of layers in radiator
     * @param n1 refractive index of the first layer (most downstream layer)
     * @param t1 thickness of the first layer
     */
    bool MakeLayers(int N, double n1, double t1);

    /**Make up focusing aerogel radiator of fixed thickness.
     * @param G sum of total radiator thickness and proximity distance
     * @param n1 refractive index of the first layer (most downstream layer)
     * @param rt1 ratio of thickness of the first layer and proximity distance
     */
    bool MakeGabarit(double G, double n1, double rt1);

    /**Make up focusing aerogel radiator with given number of layers,
     * fixed total radiator thickness, refractive index of the first layer.
     * @param N number of layers in radiator
     * @param G sum of total radiator thickness and proximity distance
     * @param n1 refractive index of the first layer (most downstream layer)
     */
    bool MakeFixed(int N, double G, double n1);

    /**Calculate radius resolution, number of photoelectrons, etc.*/
    Resolution& Calculate();

    /**Optimize focusing radiator for given number of layers, fixed total radiator thickness and refractive index of the first layer
     * @param N number of layers in radiator
     * @param G sum of total radiator thickness and proximity distance
     * @param n1 refractive index of the first layer (most downstream layer)
     */ 
    bool OptimizeNT(int N, double G, double n1);

    /**Optimize focusing radiator with a polynomial parameterization of the refractive index profile:
     *  n(x) = n(x1) + (x-x1)*P_k(x),
     * where x - coordinate along depth of the radiator starting from the face close to the PD,
     *       x1 - x value corresponding to the middle of the first layer,
     *       P_k(x) - polynomial of degree k, where k=@c npol-1. 
     * @param npol degree of polynomial, the same as the number of free coefficients during optimization. 
     * @param G sum of total radiator thickness and proximity distance
     * @param n1 refractive index of the first layer (most downstream layer)
     */ 
    bool OptimizePol(int N, int npol, double G, double nmax);

    /**Set proximity distance: distance from output face of the radiator and photodetection plane*/
    void SetProximityDistance(double d) { t0=d; optimization=0; }
    /**Set beta of particle to optimize radiator for*/
    void SetBeta(double b) { beta=(b>1||b<0)?beta=1:beta=b; optimization=0; }
    /**Set wavelength to optimize radiator for*/
    void SetWavelength(double wl) { wavelength=wl; }
    /**Set refractive index of the @c i-th layer*/
    void SetIndex(int i,double ri) { if( i<nlayers && ri>=1.0 ) { vn[i]=ri; optimization=0; } }
    /**Set thickness of the @c i-th layer*/
    void SetThickness(int i,double t) { if( i<nlayers && t>0.0 ) { vt[i]=t; optimization=0; } }
    /**Set light scattering length for aerogel*/
    void SetScatteringLength(double l) { scatteringLength=l; optimization=0; }
    /**Set photon detection efficiency spectrum*/
    void SetPDefficiency(Spectrum& eff) { pdEff = eff; optimization=0; }
    /**Set photon detector pixel size*/
    void SetPixelSize(double ps) { pixelSize = ps; optimization=0; }

    /**Get number of layers*/
    double GetNlayers() const { return nlayers; }
    /**Get proximity distance*/
    double GetProximityDistance() const { return t0; }
    /**Get beta parameter*/
    double GetBeta() const { return beta; }
    /**Get wavelength parameter*/
    double GetWavelength() const { return wavelength; }
    /**Get refractive index of the @c i-th layer, i=0..N-1*/
    double GetIndex(int i) const { return i<nlayers?vn[i]:NAN; }
    /**Get thickness of the @c i-th layer, i=0..N-1*/
    double GetThickness(int i) const { return i<nlayers?vt[i]:NAN; }
    /**Get light scattering length*/
    double GetScatteringLength() const { return scatteringLength; }
    /**Get minimum layer thickness*/
    double GetMinimumThickness() const { return *std::min_element(vt.begin(),vt.end()); }
    /**Get photon detection efficiency spectrum*/
    const Spectrum& GetPDefficiency() const { return pdEff; }
    /**Get photon detector pixel size*/
    double GetPixelSize() const { return pixelSize; }
    /**Get last optimization type*/
    int GetOptimizationType() const { return optimization; }
    /**Get last result of resolution calculation*/
    const Resolution& GetResolution() const { return result; }

    double GetAngleResolutionPerTrack() const;
    
    /**Get total thickness of radiator*/
    double GetTotalThickness() const { return std::accumulate(vt.begin(),vt.end(),0.); }

    /**Formula of refractive index expressed through the refracted Cherenkov angle tangent*/
    double rindex(double tg)
    {
        double tg2=tg*tg;
        return sqrt((1+(rn[0]*rn[0]+1)*tg2)/(tg2+1));
    }

    /**Recalculate refractive index of aerogel given at @c wl1 as @c n to @c wl2
     * Scaled from fused quartz index considering
     * (n^2-1) proportional to density*/
    static double AerogelRefIndex(double n, double wl1, double wl2);

private:
    /**Formula of tangent of the Cherenkov angle produced in the layer @c k
     * and refracted in the layer @c l*/
    double Tkl(int k, int l) {
        double nk2 = rn[k]*rn[k];
        return sqrt((nk2-1)/(rn[l]*rn[l]-nk2+1));
    }

    /**Ratio of Tkl to the tangent of the raw Cherenkov angle*/
    double Skl(int k, int l) {
        return 1/sqrt(rn[l]*rn[l]-rn[k]*rn[k]+1);
    }

    /**Compute a next layer of the focusing radiator*/
    int MakeAlayer();

    /**Transform reduced indices and thicknesses to absolute values*/
    void GoToAbs();

    /**Apply the NT parameterization to the radiator*/
    void ApplyNTParameterization(const double* p);
    /**Evaluate angle resolution with NT parameterization*/
    double EvalResolutionNT(const double* p);
    
    /**Apply the polynomial parameterization to the radiator*/
    void ApplyPolParameterization(const double* p);
    /**Evaluate angle resolution approximating density profile by polynomial*/
    double EvalResolutionPol(const double* p);

    /**Create an instance of a minimizer if it does not exist and reset minimizer if it already existed*/ 
    static void InitializeMinimizer();
};

#endif

