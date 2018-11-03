#ifndef RichGlobals_hh
# define RichGlobals_hh

# include "globals.hh"
# include "CLHEP/Units/SystemOfUnits.h"

// Conversion factor from photon wavelength in nanometers to photon energy in eV
static const G4double kPhotMomWaveConv=h_Planck*c_light/CLHEP::eV/CLHEP::nm;
//Tolerance for comparison of double values
static const G4double kTolerance = 2.2e-14;

// forward declaration of instantiated templates
const G4double& min(const G4double&,const G4double&);
const G4double& max(const G4double&,const G4double&);

#endif
