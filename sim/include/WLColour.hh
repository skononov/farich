/* Class written for Geant4 that has a constructor for monochromatic
 * light representation of given wavelength on a display specified by
 * gamma value.
 */

#ifndef WLCOLOUR_HH
# define WLCOLOUR_HH

# include "globals.hh"
# include "G4Colour.hh"
# include "ColorRep.hh"

typedef G4Colour G4Color;

class WLColour {
public:
	WLColour(G4double wl,G4double g=1.0) :
		wavelength(wl),
		gamma(g)
	{}

	G4Colour GetG4Colour(double a=1.0) {
		cr_convert(wavelength/nanometer,RGB,gamma);
		return G4Colour(RGB[0],RGB[1],RGB[2],a);
	}
	G4Colour GetG4Color(double a=1.0) {
		return GetG4Colour(a);
	}

	double GetWL() const { return wavelength; }
	double GetGamma() const { return gamma; }

	void SetWL(double wl) { wavelength=wl; }
	void SetGamma(double g) { gamma=g; }

private:
	double wavelength;
	double gamma;
	double RGB[3];
};

#endif
