#include "MLRadiatorDescription.hh"
#include <iostream>

double MLRadiatorDescription::tolerance = 1e-6;

MLRadiatorDescription::MLRadiatorDescription(double d, double b) :
	t0(d),
	beta(b),
	totalThickness(0)
{}

MLRadiatorDescription::~MLRadiatorDescription()
{
}

void MLRadiatorDescription::clear()
{
	rn.clear();
	rt.clear();
	vn.clear();
	vt.clear();
	materials.clear();
	totalThickness=0;
}

void MLRadiatorDescription::GoToAbs()
{
	int N=rn.size()-1;
	vn.resize(N);
	vt.resize(N);
	materials.resize(N);
	totalThickness=0;
	for(int i=0; i<N; i++) {
		vn[i]=rn[i+1]/beta;
		vt[i]=rt[i+1]*t0;
		materials[i]="Aerogel";
		totalThickness+=vt[i];
	}
}

void MLRadiatorDescription::AddAlayer(double ri,double t)
{
	vn.push_back(ri);
	vt.push_back(t);
	materials.push_back("Aerogel");
	totalThickness+=t;
}

void MLRadiatorDescription::AddAlayer(const char *name, double t)
{
	vn.push_back(1);
	vt.push_back(t);
	materials.push_back(name);
	totalThickness+=t;
}

int MLRadiatorDescription::MakeAlayer()
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

int MLRadiatorDescription::MakeLayers(int N, double n1, double t1)
{ //fixed number of layers. Input parameters: beta, t0, n1, t1, N.
	clear();

	rn.reserve(N+1);
	rt.reserve(N+1);

	if( N<=0 ) return 0;

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

int MLRadiatorDescription::MakeGabarit(double G, double n1, double rt1)
{ //fixed total thickness of radiator. Input parameters: beta, t0, G, n1, rt1.

	double t1 = rt1*t0; //thickness of the first layer

	double Tmax=G-t0; //maximum thickness of radiator (may be unachievable)

	clear();

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

int MLRadiatorDescription::MakeFixed(int N, double G, double n1)
{ //fixed number of layers and total thickness of radiator.
  //Input parameters: beta, t0, t0, N, G, n1
	clear();

	if( N<=0 ) return 0;

	if( N==1 ) { //single layer case
		double t1=G-t0;
		vn.push_back(n1);
		vt.push_back(t1);
		materials.push_back("Aerogel");
		totalThickness=t1;
		return 1;
	}

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

	return vn.size();
}

int MLRadiatorDescription::MakeMultiRing(int N, double  G, double n1, double n2)
{ //multi-ring option made of two series of layers of different refractive index.
  //Input parameters: t0, N, G, n1, n2
	clear();

	vn.resize(N);
	vt.resize(N);
	materials.resize(N);

	totalThickness=G-t0;
	double t1 = totalThickness/N;
	for (int i=0; i<N; i++) {
		if (i%2==0)
			vn[i]=n1;
		else
            vn[i]=n2;
		vt[i]=t1;
		materials[i]="Aerogel";
	}
	return N;
}

void MLRadiatorDescription::DistortThickness(double epst,bool even)
{
	if( vt.empty() ) return;

	totalThickness=0;

	for(size_t i=0; i<vt.size(); i++) {
		if( even )
			vt[i] *= 1+epst;
		else
			vt[i] *= 1+epst*(0.5-i%2)*2;
		totalThickness+=vt[i];
	}
}

void MLRadiatorDescription::DistortIndex(double epsd,bool even)
{
	if( vn.empty() ) return;

	for(size_t i=0; i<vn.size(); i++) {
		if( even )
			vn[i] += (vn[i]-1)*epsd;
		else
			vn[i] += (vn[i]-1)*epsd*(0.5-i%2)*2;
	}
}

int MLRadiatorDescription::DistortUniformity(double epsd,bool even)
{
    //break each layer in that amount of sublayers
	static const size_t q=5;

	if( vn.empty() ) return 0;

	int Nnew=q*vn.size();
	std::vector<double> n_new(Nnew), t_new(Nnew);
	std::vector<std::string> m_new(Nnew);

	materials.resize(Nnew);
	size_t i_new=0;

	for(size_t i=0; i<vn.size(); i++) {
		for(size_t k=0; k<q; k++) {
			if( even )
				n_new[i_new] = vn[i] + ((double)k/(q-1)-0.5)*(vn[i]-1)*epsd;
			else
				n_new[i_new] = vn[i] + ((double)k/(q-1)-0.5)*(vn[i]-1)*epsd*(0.5-i%2)*2;

			t_new[i_new] = vt[i]/q;
			m_new[i_new] = materials[i];
			i_new++;
		}
	}

	vn=n_new;
	vt=t_new;
	materials=m_new;

	return Nnew;
}

