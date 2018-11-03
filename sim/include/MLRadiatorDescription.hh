#ifndef MLRadiatorDescription_h
#define MLRadiatorDescription_h 1

# include <vector>
# include <string>
# include <cmath>

class MLRadiatorDescription {
private:
	static double tolerance; ///< allowed coordinate shift calculation error

	double t0;         ///< proximity distance
	double beta;       ///< particle velocity that we need to optimize for

	std::vector<std::string> materials; ///< material name of layers
	std::vector<double> vn;  ///< indices of refraction
	std::vector<double> vt;  ///< thicknesses in usual units

	std::vector<double> rn;  ///< reduced indices for calculations
	std::vector<double> rt;  ///< reduced thicknesses for calculations

	double totalThickness; ///< total thickness of layers in usual units

	double T10; ///< auxilary coefficient

public:
	/**Set tolerance - calculation accuracy parameter*/
	static void SetTolerance(double eps) { tolerance=eps; }

	/**Constructor.
	 * @param d proximity distance
	 * @param b beta parameter
	 */
	MLRadiatorDescription(double d=1.0, double b=1.0);
	/**Destructor*/
	~MLRadiatorDescription();

	/**Add a layer of aerogel with refractive index @c ri and thickness @c t.*/
	void AddAlayer(double ri,double t);
	/**Add a layer of material @c name with thickness @c t*/
	void AddAlayer(const char* name,double t);

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

	/**Make up "multi-ring" aerogel radiator of layers with two series of refractive index.
	 * @param N number of layers in radiator
	 * @param G sum of total radiator thickness and proximity distance
	 * @param n1 index of the first series
	 * @param n2 index of the second series
	 * Thickness of the radiator is determined as difference between @c G and proximity distance.
	 */
	int MakeMultiRing(int N, double  G, double n1, double n2);

	/**Distort thicknesses of every layer by given relative magnitude @c epst.
	 * If @c even is true, distort every layer in the same direction, otherwise distort
	 * thickness of every other layer in opposite direction.
	 * For studying stability of focusing condition in real aerogel production.
	 */
	void DistortThickness(double epst, bool even);

	/**Distort refractive indices of every layer by given relative density magnitude @c epsd.
	 * If @c even is true, distort every layer in the same direction, otherwise distort
	 * index of every other layer in opposite direction.
	 */
	void DistortIndex(double epsd, bool even);

	/**Introduce linear refractive index non-uniformity in each layer of @c epsd
	 * relative density magnitude. It is done by means of breaking each layers into several
	 * sub-layers with indices varying linearly. Function returns number of layers in resulting
	 * radiator. If @c even is true, index gradient in every layer goes in the same direction,
	 * otherwise every other layer has opposite gradient direction.
	 */
	int DistortUniformity(double epsd, bool even);

	/**Reset radiator*/
	void clear();

	/**Set proximity distance: distance from output face of the radiator and photodetection plane*/
	void SetProximityDistance(double d) { t0=d; }
	/**Set beta of particle to calculate radiator for*/
	void SetBeta(double b) { beta=(b>1||b<0)?1:b; }

	/**Get proximity distance*/
	double GetProximityDistance() const { return t0; }
	/**Get beta parameter*/
	double GetBeta() const { return beta; }

	/**Get total thickness of radiator*/
	double GetTotalThickness() const { return totalThickness; }
	/**Get number of layers in radiator*/
	int GetNlayers() const { return vn.size(); }
	/**Get refractive index of the @c i-th layer, i=0..N-1*/
	double GetIndex(int i) const { return vn[i]; }
	/**Get thickness of the @c i-th layer, i=0..N-1*/
	double GetThickness(int i) const { return vt[i]; }
	/**Get material name of the @c i-th layer, i=0..N-1*/
	const char* GetMaterialName(int i) const { return materials[i].c_str(); }

	/**Set refractive index of the @c i-th layer*/
	void SetIndex(int i,double ri) { vn[i]=ri; }

private:
	/**Compute a next layer of the focusing radiator*/
	int MakeAlayer();

	/**Transform reduced indices and thicknesses to absolute values*/
	void GoToAbs();

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
};

#endif

