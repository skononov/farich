#ifndef RichRootEvent_hh
# define RichRootEvent_hh

# include <cmath>
# include "TObject.h"
# include "TString.h"
# include "TClonesArray.h"


class RichRootParticle //: public TObject
{
public:
    TString name;
	Float_t x0, y0, z0; //point of origin [mm]
	Float_t vx, vy, vz; //direction at origin [unit vector]
    Float_t p; //momentum [GeV/c]
    Float_t m; //mass [GeV/c^2]

	RichRootParticle() : //TObject(),
    	name("undef"),
		x0(0.), y0(0.), z0(0.), vx(0.), vy(0.), vz(0.),
    	p(0), m(0)
	{}

	~RichRootParticle() {}

	Double_t Beta() const { return p==0.?0.:1/sqrt(1+m*m/(p*p)); }

	ClassDef(RichRootParticle,1)
};

class RichRootHit : public TObject
{
public:
	Int_t   pxid; //hit pixel ID
	Float_t x, y, z; //hit position by center of pixel [mm]
	Float_t xe, ye, ze; //hit exact position [mm]
	Float_t x0, y0, z0; //origin of the photon [mm]
	Float_t vx, vy, vz; //photon direction at origin w.r.t. lab frame
	Float_t wl; //photon wavelength [nm]
	Float_t thetac0, phic0; //MC truth photon angles, [radian]
	Float_t radius; //reconstructed radius [mm]
	Float_t thetac, phic; //reconstructed photon angles [radian]
	Bool_t scattered; //if the photon was previously scattered
	Short_t n; //number of photons detected by the pixel
	Short_t layer; //layer of photon origin

	RichRootHit() : //TObject(),
		pxid(0), x(0.), y(0.), z(0.), xe(0.), ye(0.), ze(0.),
		vx(0.), vy(0.), vz(0.), wl(0.), thetac0(0.), phic0(0.),
		radius(0.), thetac(0.), phic(0.),
        scattered(0), n(0), layer(0)
	{}

	~RichRootHit() {}

	ClassDefNV(RichRootHit,1)
};

class RichRootEvent : public TObject
{
private:
    UInt_t capacity;

public:
	UInt_t nhits;
	RichRootParticle particle; //particle
	TClonesArray *hits; //-> array of hits

	static TClonesArray *ghits; //-> array of hits

	RichRootEvent();
	~RichRootEvent() { hits->Clear(); }

	RichRootHit* AddHit();

	void Clear(const Option_t* opt="") { hits->Clear(); nhits=0; }

	static void Reset() { if( ghits ) { ghits->Delete(); delete ghits; } }

	ClassDefNV(RichRootEvent,1)
};

#endif

