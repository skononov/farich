#ifndef RichPMT_h
# define RichPMT_h 1
# include <set>

# include "globals.hh"
# include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class RichHitsCollection;

class RichPMT : public G4VSensitiveDetector {
public:
	RichPMT(G4String DetName);
	virtual ~RichPMT();

	void SetPixelSizeAndSpacing(G4double size,G4double spacing,G4bool round=false);

	void Initialize(G4HCofThisEvent *);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *);
	void EndOfEvent(G4HCofThisEvent *);

	G4bool ProcessHits_constStep(const G4Step *aStep);

	void clear();
	void DrawAll();
	void PrintAll();

private:
	G4int GetPixelAndRound(G4double &x, G4double &y);

private:
	RichHitsCollection *collection;
	G4int collectionID;
	G4int NumHits;

	G4double pixelSize;
	G4double pixelSpacing;
	G4bool   roundPixel;
	G4double pixelSpacingX;
	G4double pixelSpacingY;

	G4double x0, y0; //position of the first pixel left lower corner
	G4int nxPixels, nyPixels; //number of pixels along X & Y

	std::set<G4int> pixelCollection;
};

#endif
