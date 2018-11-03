#ifndef RichHitInfo_hh
# define RichHitInfo_hh

struct HitInfo_0 {
	UInt_t event; //event number in this run
	UInt_t layer; //layer of origin
	Float_t x, y; //hit position by pixel, mm
	Float_t x0, y0, z0; //vertex position, mm
	Float_t theta0, phi0; //initial angles of photon direction in point relative to particle direction, degrees
	Float_t theta, phi; //angles of photon direction at photodetector relative to detector's normal, degrees
	Float_t wl; //wavelength of the photon in nm
	Float_t radius, angle, beta; //reconstructed values
	Float_t xe, ye; //hit exact position, mm
	UInt_t npe; //number of photoelectrons in the pixel
	Bool_t scattered; //if the photon was previously scattered
};

static const char * const hitinfofmt_0="event/i:layer:x/F:y:x0:y0:z0:theta0:phi0:theta:phi:wl:radius:angle:beta:xe:ye:npe/i:scattered/O";

struct HitInfo {
	UInt_t event; //event number in this run
	UInt_t layer; //layer of origin
	Float_t x, y; //hit position by pixel, mm
	Float_t x0, y0, z0; //vertex position, mm
	Float_t theta0, phi0; //actual angles of photon emission direction in the particle's frame, degrees
	Float_t theta, phi; //actual angles of photon direction at photodetector plane relative to lab frame, degrees
	Float_t wl; //wavelength of the photon in nm
	Float_t radius; //reconstructed ring radius, mm
	Float_t angle, alpha; //reconstructed polar and azimutal angle of photon: angle [mrad], alpha [degrees]
	Float_t beta; //reconstructed value of particle's beta
	Float_t xe, ye; //hit exact position, mm
	UInt_t npe; //number of photoelectrons in the pixel
	Bool_t scattered; //if the photon was previously scattered
};

static const char * const hitinfofmt="event/i:layer:x/F:y:x0:y0:z0:theta0:phi0:theta:phi:wl:radius:angle:alpha:beta:xe:ye:npe/i:scattered/O";
#endif
