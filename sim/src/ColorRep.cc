/**************************************************************************
*           RGB representation of visible wavelengths.                    *
* Conversion from the original fortran code of Dan Bruton (astro@tamu.edu)*
* The original can be found at                                            *
* http://www.physics.sfasu.edu/astro/color.html                           *
* Created by Sergey Kononov, 2004                                         *
***************************************************************************/
#include <iostream>
#include <cmath>

//Function converts the given wavelength in nanometers
//to RGB relative intensities for a specific gamma.
int cr_convert(double wl,double *rgb,double gamma=1.0)
{
	double r, g, b;
	double si;

	if (gamma<0.) {
		std::cerr << "Illegal value of gamma " << gamma << std::endl;
		return -2;
	}
	if (wl<0.) {
		std::cerr << "Illegal value of wavelength " << wl << std::endl;
		return -3;
	}

	if (wl<380) {
		r = 1.;
		g = 0.;
		b = 1.;
	}
	else if (wl<=440) {
		r = (440.-wl)/(440.-380.);
		g = 0.;
		b = 1.;
	}
	else if (wl<=490) {
		r = 0.;
		g = (wl-440.)/(490.-440.);
		b = 1.;
	}
	else if (wl<=510) {
		r = 0.;
		g = 1.;
		b = (510.-wl)/(510.-490.);
	}
	else if (wl<=580) {
		r = (wl-510.)/(580.-510.);
		g = 1.;
		b = 0.;
	}
	else if (wl<=645) {
		r = 1.;
		g = (645.-wl)/(645.-580.);
		b = 0.;
	}
	else {
		r = 1.;
		g = 0.;
		b = 0.;
	}

// LET THE INTENSITY SI FALL OFF NEAR THE VISION LIMITS
	if (wl>700.) 		si = .3+.7*(780.-wl)/(780.-700.);
	else if (wl<420.) 	si = .3+.7*(wl-380.)/(420.-380.);
	else                si = 1.;

	if (si<0.3) {
		si = 0.3;
/*		if (wl<380.) //UltraViolet
			return -1;
		else         //InfraRed
			return 1;
*/	}

// GAMMA ADJUST AND WRITE IMAGE TO AN ARRAY
	rgb[0] = pow(si*r,gamma);
	rgb[1] = pow(si*g,gamma);
	rgb[2] = pow(si*b,gamma);

	return 0;
}

