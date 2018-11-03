/**************************************************************************
*           RGB representation of visible wavelengths.                    *
* Conversion from the original fortran code of Dan Bruton (astro@tamu.edu)*
* The original can be found at                                            *
* http://www.physics.sfasu.edu/astro/color.html                           *
* Created by Sergey Kononov, 2004                                         *
* Header file. Declares cr_convert function converting wavelength to RGB  *
***************************************************************************/
#ifndef COLORREP_HH
# define COLORREP_HH

extern int cr_convert(double wl, double *rgb, double gamma=1.0);

#endif

