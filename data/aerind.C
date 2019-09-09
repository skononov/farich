#include "TF1.h"

double AerogelRefIndex(double* x, double *par)
{
    // Parameters of the one-pole Sellmeier formula fit to Novosibirsk aerogel n=1.03 data.
    // T. Bellunato et al., "Refractive index dispersion law of silica aerogel",
    // Eur. Phys. J. C 52 (2007) 759-764

	double wl = *x, n0 = par[0], wl0 = par[1];
    const double LHCb_a0 = 0.05639, LHCb_wl0sqr = 83.22 * 83.22;
    double LHCb_RI2m1ref = LHCb_a0 / (1 - LHCb_wl0sqr / (wl0 * wl0)); //(n**2-1) of LHCb aerogel at wl0
    double ri2m1_lhcb = LHCb_a0 / (1 - LHCb_wl0sqr / (wl * wl));      //(n**2-1) of LHCb aerogel at wl

    return sqrt(1 + (n0 * n0 - 1) / LHCb_RI2m1ref * ri2m1_lhcb);
}

TF1* faerind = new TF1("faerind",AerogelRefIndex,250.,700.,2);

