#include "RichPMT.hh"

#include <cmath>

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "Randomize.hh"

#include "RichGlobals.hh"
#include "RichParameters.hh"
#include "RichSteppingAction.hh"
#include "RichHit.hh"
#include "RichUserTrackInformation.hh"

#include "G4SystemOfUnits.hh"

using std::pair;
using std::set;

RichPMT::RichPMT(G4String DetName) :
	G4VSensitiveDetector(DetName),
	collectionID(-1),
	NumHits(0)
{
    G4String HCname;
	collectionName.insert(HCname="RichHitsCollection");
}

RichPMT::~RichPMT() { }

void RichPMT::Initialize(G4HCofThisEvent* HCE)
{
	collection =
		new RichHitsCollection(SensitiveDetectorName,collectionName[0]);
	if (collectionID < 0) {
		collectionID =
			G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	}
	HCE->AddHitsCollection(collectionID,collection);
	NumHits=0;

	pixelCollection.clear();
}

G4bool RichPMT::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
	if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhoton())
		G4cerr << "RichPMT::ProcessHits() called directly!!!\a\n";
	return ProcessHits_constStep(aStep);
}

void RichPMT::EndOfEvent(G4HCofThisEvent *) {}

void RichPMT::clear() {
	pixelCollection.clear();
}

void RichPMT::DrawAll() {}

void RichPMT::PrintAll() {}

void RichPMT::SetPixelSizeAndSpacing(G4double size,G4double spacing,G4bool round)
{
	pixelSize=size;
	pixelSpacing=spacing;
	roundPixel=round;

	if( round ) {
/*		pixelSpacingX=pixelSpacing;
		pixelSpacingY=0.5*sqrt(3)*pixelSpacing;
		nxPixels=int((2*kDetHalfX-pixelSize-0.5*pixelSpacingX)/pixelSpacingX)+1;
		nyPixels=int((2*kDetHalfY-pixelSize)/pixelSpacingY)+1;
        x0=
*/	} else {
		nxPixels=int(2*kDetHalfX/pixelSpacing);
		nyPixels=int(2*kDetHalfY/pixelSpacing);
		x0=-nxPixels/2.*pixelSpacing;
		y0=-nyPixels/2.*pixelSpacing;
	}
}

G4int RichPMT::GetPixelAndRound(G4double &x, G4double &y)
{
	if( pixelSpacing<=0 ) return 0;

	if( roundPixel ) {
// TODO
	} else {
		G4double xr=fmod(x-x0,pixelSpacing)-pixelSpacing/2;
		G4double yr=fmod(y-y0,pixelSpacing)-pixelSpacing/2;

		if( fabs(xr)>pixelSize/2. || fabs(yr)>pixelSize/2. ) return -1; //not hit

		G4int idx = int((x-x0)/pixelSpacing);
		G4int idy = int((y-y0)/pixelSpacing);
		G4int id  = idy*nxPixels+idx;

		x=x0+(idx+0.5)*pixelSpacing;
		y=y0+(idy+0.5)*pixelSpacing;

		return id;
	}
	return 0;
}


G4bool RichPMT::ProcessHits_constStep(const G4Step *aStep)
{
	G4Track *aTrack = aStep->GetTrack();

	if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhoton())
		return false;

	RichUserTrackInformation *trackInformation =
		(RichUserTrackInformation *) aTrack->GetUserInformation();
	trackInformation->AddTrackStatusFlag(trackInformation->hitPMT);

	G4int pixelId=0;
	G4ThreeVector expos=aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector pos=expos;
	if( pixelSpacing>0. ) {
		pixelId=GetPixelAndRound(pos[0],pos[1]);

		if( pixelId<0 ) return false; //not hit

		pair<set<G4int>::iterator, bool> res=pixelCollection.insert(pixelId);

		if( !res.second ) { // this pixel has already been hit
			for(size_t i=0; i<collection->GetSize(); i++) {
				RichHit* aHit=(RichHit*)collection->GetHit(i);
				if( aHit->GetPixelId()==pixelId ) {
					aHit->Hit(); //increment hit count for the pixel
					if (verboseLevel>1)
						G4cout << aHit->GetNhits() << " photons in the pixel " << pixelId << G4endl;
					break;
				}
			}
			return true;
		}
	}

	G4double wavelength = kPhotMomWaveConv/
		(aTrack->GetTotalEnergy()/eV)*nanometer;

	RichHit *aHit = new RichHit();
	aHit->SetPixelId(pixelId);
	aHit->SetPosition(pos);
	aHit->SetExactPosition(expos);
	aHit->SetVertexPos(aTrack->GetVertexPosition());
	aHit->SetDirection(aStep->GetPreStepPoint()->GetMomentumDirection());
	aHit->SetVertexDirection(aTrack->GetVertexMomentumDirection());
	aHit->SetWL(wavelength);
	aHit->SetLayerNo(trackInformation->GetLayerNo());

	G4StepPoint *thePostPoint = aStep->GetPostStepPoint();
	G4String pname = thePostPoint->GetProcessDefinedStep()->GetProcessName();

	if (trackInformation->GetScatterCount()==0 &&
		trackInformation->GetReflectionCount()==0 &&
		trackInformation->IsFromPrimary() ) {
		if (verboseLevel>1)
			G4cout << "Clean hit" << G4endl;
		aHit->Hit();
	}

	collection->insert(aHit);
	NumHits++;

	return aHit->IsHit();
}
