#ifndef RichParticlePropagation_h
# define RichParticlePropagation_h 1

# include "Rtypes.h"
# include "globals.hh"
# include "G4ThreeVector.hh"

class TF1;

class RichParticlePropagation
{
public:
	RichParticlePropagation(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4double);
	~RichParticlePropagation() {}

	void Initialize(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4double);

	G4bool Propagate(G4double, G4ThreeVector &, G4ThreeVector &);

	void Draw() const;

private:
	G4bool inited, straight;
    TF1 *funcZshift;
	G4double radius, vPerp, vPar;
    G4double hmax, hperiod;
    G4ThreeVector position, tangent, perpTangent, fieldDir, curvature;

    static const Double_t approx_limit;
    static const Double_t delta;

	static Double_t EvalZshift(Double_t *,Double_t *);
};

#endif
