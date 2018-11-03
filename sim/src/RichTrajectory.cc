#include <G4ios.hh>
#include <G4TrajectoryPoint.hh>
#include <G4Trajectory.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4ThreeVector.hh>
#include <G4Polyline.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4VVisManager.hh>
#include <G4Polymarker.hh>
#include <G4ParticleDefinition.hh>

#include "RichTrajectory.hh"
#include "RichGlobals.hh"
#include "WLColour.hh"

#include "G4SystemOfUnits.hh"

G4Allocator < RichTrajectory > RichTrajectoryAllocator;

RichTrajectory::RichTrajectory() :
	G4Trajectory(), drawit(false), wavelength(0.0)
{
	particleDefinition = 0;
}

RichTrajectory::RichTrajectory(const G4Track * aTrack) :
	G4Trajectory(aTrack), drawit(false), wavelength(0.0)
{
	particleDefinition = aTrack->GetDefinition();
    if (particleDefinition == G4OpticalPhoton::OpticalPhoton())
		wavelength = kPhotMomWaveConv/(aTrack->GetTotalEnergy()/eV)*nanometer;
}

RichTrajectory::RichTrajectory(RichTrajectory & right) :
	G4Trajectory(right), drawit(right.drawit)
{
	particleDefinition = right.particleDefinition;
	wavelength = right.wavelength;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichTrajectory::~RichTrajectory() {}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichTrajectory::DrawTrajectory(G4int i_mode) const
{
	//Taken from G4VTrajectory and modified to select colours based on particle
	//type and to selectively eliminate drawing of certain trajectories.

	if (!drawit) return;

	// If i_mode>=0, draws a trajectory as a polyline and, if i_mode!=0,
	// adds markers - yellow circles for step points and magenta squares
	// for auxiliary points, if any - whose screen size in pixels is
	// given by abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
	// visible markers.

	G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
	if (!pVVisManager) return;

	const G4double markerSize = abs(i_mode) / 1000;
	G4bool lineRequired(i_mode >= 0);
	G4bool markersRequired(markerSize > 0.);

	G4Polyline trajectoryLine;
	G4Polymarker stepPoints;
    G4Polymarker auxiliaryPoints;

	for (G4int i = 0; i < GetPointEntries(); i++) {
		G4VTrajectoryPoint *aTrajectoryPoint = GetPoint(i);

		const std::vector<G4ThreeVector>* auxiliaries
			= aTrajectoryPoint->GetAuxiliaryPoints();
		if (auxiliaries) {
			for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
				const G4ThreeVector pos((*auxiliaries)[iAux]);
				if (lineRequired)
					trajectoryLine.push_back(pos);
				if (markersRequired)
					auxiliaryPoints.push_back(pos);
			}
		}

		const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
		if (lineRequired)
			trajectoryLine.push_back(pos);
		if (markersRequired)
			stepPoints.push_back(pos);
	}

	if (lineRequired) {
		G4Colour colour;

/*		if (wavelength!=0) { //if photon
			WLColour wlcol(wavelength); //draw it of corresponding colour
			colour = wlcol.GetG4Color();
		}
		else*/ if (particleDefinition->GetPDGCharge()>0)
			colour = G4Colour(1., 0., 0.); //if positive draw red
		else if (particleDefinition->GetPDGCharge()<0)
			colour = G4Colour(0., 0., 1.); //if negative draw blue
		else if (particleDefinition->GetPDGCharge()==0)
			colour = G4Colour(1., 1., 0.); //if neutral draw yellow

		G4VisAttributes trajectoryLineAttribs(colour);
		trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
		pVVisManager->Draw(trajectoryLine);
	}

	if (markersRequired) {
		auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
		auxiliaryPoints.SetScreenSize(markerSize);
		auxiliaryPoints.SetFillStyle(G4VMarker::filled);
		G4VisAttributes auxiliaryPointsAttribs(G4Colour(.5,.5,.5));  //Grey
		auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
		pVVisManager->Draw(auxiliaryPoints);

		stepPoints.SetMarkerType(G4Polymarker::circles);
		stepPoints.SetScreenSize(markerSize);
		stepPoints.SetFillStyle(G4VMarker::filled);
		G4VisAttributes stepPointsAttribs(G4Colour(1., 1., 1.));	//White
		stepPoints.SetVisAttributes(&stepPointsAttribs);
		pVVisManager->Draw(stepPoints);
	}
}
