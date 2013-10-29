#ifndef MLADescription_h
#define MLADescription_h 1

# include <vector>
# include <string>
# include <cmath>

/**Derived from MLRadiatorDescription class in FARICH simulation.
   Date: 27.02.2009
*/

class Spectrum;

static const int Nr=100;

struct MLAResult {
	bool valid;
	double npe;
	double radius;
	double sigma1;
	double sigma_t;
	double rmin, rmax;
	double r[Nr];
	double s[Nr];
};

class MLADescription {
private:
	static double tolerance; ///< allowed coordinate shift calculation error

	int nlayers;              ///< number of layers in radiator
	double t0;                ///< proximity distance
	double beta;              ///< particle velocity that we need to optimize for
	double wavelength;        ///< wavelength at which we need to optimize for, nm
	double scatteringLength;  ///< light scattering length in usual units at 400 nm

	std::vector<double> vn;   ///< indices of refraction at 400 nm
	std::vector<double> vt;   ///< thicknesses in usual units

	std::vector<double> rn;  ///< reduced indices for calculations
	std::vector<double> rt;  ///< reduced thicknesses for calculations

	double totalThickness;   ///< total thickness of layers in usual units

	double T10;              ///< auxilary coefficient

	MLAResult result;        ///< last result of calculation

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

	/**Add a layer of aerogel with refractive index @c ri, thickness @c t*/
	void AddAlayer(double ri,double t);

	/**Make up focusing aerogel radiator with given number of layers.
	 * @param N number of layers in radiator
	 * @param n1 refractive index of the first layer (most downstream layer)
	 * @param t1 thickness of the first layer
	 */
	int MakeLayers(int N, double n1, double t1);

	/**Make up focusing aerogel radiator of fixed thickness.
	 * @param G sum of total radiator thickness and proximity distance
	 * @param n1 refractive index of the first layer (most downstream layer)
	 * @param rt1 ratio of thickness of the first layer and proximity distance
	 */
	int MakeGabarit(double G, double n1, double rt1);

	/**Make up focusing aerogel radiator with given number of layers and
	 * fixed total radiator thickness.
	 * @param N number of layers in radiator
	 * @param G sum of total radiator thickness and proximity distance
	 * @param n1 refractive index of the first layer (most downstream layer)
	 */
	int MakeFixed(int N, double G, double n1);

	/**Reset radiator*/
	void clear();

	/**Print radiator configuration*/
	void Print(const char* pfx="") const;

	/**Calculate radius resolution, number of photoelectrons, etc. of the radiator
	 * for given photon detector efficiency @c eff, pixel size @c ps, and particle velocity @c b*/
	MLAResult& Calculate(Spectrum& eff,double ps=0,double b=0);

	/**Set proximity distance: distance from output face of the radiator and photodetection plane*/
	void SetProximityDistance(double d) { t0=d; }
	/**Set beta of particle to optimize radiator for*/
	void SetBeta(double b) { beta=(b>1||b<0)?beta=1:beta=b; }
	/**Set wavelength to optimize radiator for*/
	void SetWavelength(double wl) { wavelength=wl; }
	/**Set refractive index of the @c i-th layer*/
	void SetIndex(int i,double ri) { if( i<nlayers && ri>=1.0 ) vn[i]=ri; }
	/**Set thickness of the @c i-th layer*/
	void SetThickness(int i,double t) { if( i<nlayers && t>0.0 ) { vt[i]=t; } }
	/**Set light scattering length for all the layers*/
	void SetScatteringLength(double l) { scatteringLength=l; }
	/**Recalculate total thickness of radiator*/
	double CalculateTotalThickness();

	/**Get number of layers*/
	double GetNlayers() const { return nlayers; }
	/**Get proximity distance*/
	double GetProximityDistance() const { return t0; }
	/**Get beta parameter*/
	double GetBeta() const { return beta; }
	/**Get wavelength parameter*/
	double GetWavelength() const { return wavelength; }
	/**Get total thickness of radiator*/
	double GetTotalThickness() const { return totalThickness; }
	/**Get refractive index of the @c i-th layer, i=0..N-1*/
	double GetIndex(int i) const { return i<nlayers?vn[i]:NAN; }
	/**Get thickness of the @c i-th layer, i=0..N-1*/
	double GetThickness(int i) const { return i<nlayers?vt[i]:NAN; }
	/**Get light scattering length of the @c i-th layer, i=0..N-1*/
	double GetScatteringLength() const { return scatteringLength; }
	/**Get minimum layer thickness*/
	double GetMinimumThickness() const;

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
	/**Compute a next layer of the focusing radiator*/
	int MakeAlayer();

	/**Transform reduced indices and thicknesses to absolute values*/
	void GoToAbs();
};

#endif

