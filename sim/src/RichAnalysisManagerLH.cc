#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>

#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4HCofThisEvent.hh"

//ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "RichGlobals.hh"
#include "RichHitInfo.hh"
#include "RichHit.hh"
#include "RichAnalysisMessenger.hh"
#include "RichDetectorConstruction.hh"
#include "RichPrimaryGeneratorAction.hh"
#include "RichAnalysisManagerLH.hh"

using namespace std;

#define SUM2_to_RMS(sum2,mean,n) (sum2<n*mean*mean?0:sqrt(sum2/n-mean*mean))


static Int_t nwLayer=0;

RichAnalysisManagerLH::RichAnalysisManagerLH() : RichVAnalysisManager()
{
	nRadii=0;
	vRadii=0;
	vLayerRadii=0;
}

RichAnalysisManagerLH::~RichAnalysisManagerLH()
{}

void RichAnalysisManagerLH::book(const G4Run* run)
{
	RichVAnalysisManager::book(run);

//Initialize sum values per run
	runNumPe=0;

	runCkvAngle=0;
	run1peAngleRMS=0;
	runTrackAngleRMS=0;
	runRingRadius=0;
	run1peRadiusRMS=0;
	runTrackRadiusRMS=0;
	runBeta=0;
	run1peBeta=0;
	run1peBetaRMS=0;
	runTrackBetaRMS=0;

	runNumEvWithHit=0;
	runNumHits=0;

	nwLayer=0;

//Allocate space for values per event
	nRadii=(Int_t)(5*expectedNpe);
	vRadii=new Double_t[nRadii];
	vLayerRadii=new Double_t*[nLayers];
	for(Int_t l=0; l<nLayers; l++)
		vLayerRadii[l]=new Double_t[nRadii];
}

#include "RichAnalysisModuleLH.C"

void RichAnalysisManagerLH::finish()
{
//Free space after run
	if (layerNumPe) delete [] layerNumPe;
	if (vRadii) { delete [] vRadii; nRadii=0; }
	if (vLayerRadii) {
		for(Int_t l=0; l<nLayers; l++)
			delete [] vLayerRadii[l];
		delete [] vLayerRadii;
	}

// Some results
	runNumPe /= currentRun->GetNumberOfEvent();

	runCkvAngle      /= runNumEvWithHit;
	runRingRadius    /= runNumEvWithHit;
	run1peAngleRMS    = SUM2_to_RMS(run1peAngleRMS,runCkvAngle,runNumHits);
	runTrackAngleRMS  = SUM2_to_RMS(runTrackAngleRMS,runCkvAngle,runNumEvWithHit);
	run1peRadiusRMS   = SUM2_to_RMS(run1peRadiusRMS,runRingRadius,runNumHits);
	runTrackRadiusRMS = SUM2_to_RMS(runTrackRadiusRMS,runRingRadius,runNumEvWithHit);
	runBeta         /= runNumEvWithHit;
	run1peBeta      /= runNumHits;
	run1peBetaRMS    = SUM2_to_RMS(run1peBetaRMS,run1peBeta,runNumHits);
	runTrackBetaRMS  = SUM2_to_RMS(runTrackBetaRMS,runBeta,runNumEvWithHit);

	Int_t nwLayers=0;
	for (G4int i=0; i<nLayers; i++)
		if (hLayer1peRadius[i]->GetEntries()>0) nwLayers++;

	ancout << "------------------ Run " << currentRun->GetRunID() << " summary --------------------------" << "\n"
		   << "Number of events with at list one hit in this run " << runNumEvWithHit << "\n"
	 	   << " Number of working layers: " << nwLayers << "\n"
	 	   << " Mean number of photoelectrons: " << runNumPe << "\n"
		   << " Cherenkov angle, mrad: " << runCkvAngle/mrad << " measured, (1st layer"
		   << layerAngles[0]/mrad << " est)\n"
		   << " Cherenkov ring radius, mm: " << runRingRadius/mm << "mm\n"
		   << " Angular resolution, mrad: " << run1peAngleRMS/mrad << " per photon, "
		   << runTrackAngleRMS/mrad << " per track\n"
		   << " Radius resolution, mm: " << run1peRadiusRMS/mm << " per photon, "
		   << runTrackRadiusRMS/mm << " per track\n"
		   << " Velocity: " << runBeta << " (" << beta << "), diff=" << fabs(runBeta-beta) << "\n"
		   << " Velocity resolution: " << run1peBetaRMS << " per photon, "
		   << runTrackBetaRMS << " per track\n"
		   << "------------------------------------------------------------" << G4endl;

	ancout<<"Hits with wrongly defined layers "<<100.*nwLayer/runNumHits<<endl;

	RichVAnalysisManager::finish();
}

void RichAnalysisManagerLH::LoadBetaData()
{
	betaDataExist=false;

	if (betaDataFileName.empty()) return;

	Int_t rc=loadcal(betaDataFileName);

	if (!rc) {
		G4cerr << "Calibration data incorrect and will not be used" << G4endl;
		return;
	}
	if (Nr!=nRings) {
		G4cerr << "Current number of rings does not match the number read from file" << G4endl;
		return;
	}

	betaDataExist=true;
}


void RichAnalysisManagerLH::EndOfEventAnalysis(const RichHitsCollection *RHC)
{
	iTimer->Stop();

	if (verboseLevel>1) {
		ancout << *iTimer << "\n"
		       << " Number of Cherenkov photons " << NumCkvPhot << "\n"
		       << " Number of photons at PMT " << NumPhotAtPMT << G4endl;
	}

	hNumCkvPhot->Fill(NumCkvPhot);
	hNumPhotAtPMT->Fill(NumPhotAtPMT);

	if (RHC==0) return;

	G4double trackCkvAngle=0, trackRadius=0;
	G4int nhits=RHC->entries();

	if (!nhits) return;

	if( nhits>nRadii ) {
		delete [] vRadii;
		nRadii=2*nhits;
		vRadii=new Double_t[nRadii];
		for(Int_t l=0; l<nLayers; l++) {
			delete [] vLayerRadii[l];
			vLayerRadii[l]=new Double_t[nRadii];
		}
	}

	memset(layerNumPe,0,nLayers*sizeof(G4int));

	for (G4int ih=0; ih<nhits; ih++ )
	{
		RichHit* aHit = (RichHit*)RHC->GetHit(ih);

		if (!aHit->IsHit()) continue;

		NumPhotoelectrons++;

		G4ThreeVector hitPos = aHit->GetPosition();
		G4double x=hitPos.x()/mm, y=hitPos.y()/mm;

		hXYhits->Fill(x,y);

		G4double wl = aHit->GetWL()/nanometer;
		hPeWL->Fill(wl);

		G4int iLayer = aHit->GetLayerNo();

		//Cherenkov angle fine estimation. Account for refraction.
		G4double thetac, alpha;
		find_thetac(hitPos,iLayer,thetac,alpha);
		if( TMath::IsNaN(thetac) ) thetac=0;

		//Ring radius calculation
		G4double radius = (hitPos-imageCenter).mag();

		G4double thetac_mr = thetac/mrad, radius_mm=radius/mm;

		trackCkvAngle += thetac;
		trackRadius += radius;

		run1peAngleRMS += thetac*thetac;
		run1peRadiusRMS += radius*radius;

		hCkvAngle->Fill(thetac_mr);
		hAlphaAngle->Fill(alpha/deg);
		h1peRadius->Fill(radius_mm);

		hLayerCkvAngle[iLayer]->Fill(thetac_mr);
		hLayer1peRadius[iLayer]->Fill(radius_mm);
		layerNumPe[iLayer]++;

		Double_t b=0;

		if (betaDataExist) {
			//Beta calculation
			Int_t defLayer=iLayer;
			b = findOptimalBeta(radius,beta,defLayer);
			if( defLayer!=iLayer )
				nwLayer++;
			run1peBeta+=b;
			run1peBetaRMS+=b*b;
			h1peBeta->Fill(b);
			hLayer1peBeta[iLayer]->Fill(b);
			vRadii[NumPhotoelectrons-1] = radius;
			vLayerRadii[iLayer][layerNumPe[iLayer]-1] = radius;
		}

		//Copy hit information to the storable structure
		hitinfo.event=aHit->GetEvent();
		hitinfo.x=x;
		hitinfo.y=y;
		hitinfo.layer=iLayer+1;

		G4ThreeVector vertexPos(aHit->GetVertexPos());
		hitinfo.x0=vertexPos.x()/mm;
		hitinfo.y0=vertexPos.y()/mm;
		hitinfo.z0=vertexPos.z()/mm;

		G4ThreeVector dir0(aHit->GetVertexDirection());
		//transform dir0 to cf bound to primary particle
		dir0.rotateZ(-primDirection.phi());
		dir0.rotateY(-primDirection.theta());
		hitinfo.theta0=dir0.theta()/deg;
		hitinfo.phi0=dir0.phi()/deg;

		G4ThreeVector dir(aHit->GetDirection());
		hitinfo.theta=dir.angle(primDirection)/deg;
		hitinfo.phi=dir.phi()/deg;

		hitinfo.wl=wl;

		hitinfo.radius=radius_mm;
		hitinfo.angle=thetac_mr;
		hitinfo.alpha=alpha/deg;
		hitinfo.beta=b;

		//Fill tree for this fit
		tHits->Fill();

        //Print hit info
		if (verboseLevel>2) {
			ancout << "Hit#" << NumPhotoelectrons << "\n"
				   << " hitPos = " << hitPos/mm << "\n"
				   << " wavelength = " << wl << "nm\n"
				   << " layer# = " << iLayer+1 << "\n"
				   << " radius = " << radius_mm << "mm\n"
				   << " thetac = " << thetac_mr << "mrad\n"
				   << " beta = " << b << G4endl;
		}
	}

	Double_t trackBeta=0;

	if (NumPhotoelectrons>0) {
		runNumHits += NumPhotoelectrons;
		runNumEvWithHit++;

		runNumPe += NumPhotoelectrons;
		trackCkvAngle /= NumPhotoelectrons;
		trackRadius /= NumPhotoelectrons;

		runCkvAngle += trackCkvAngle;
		runTrackAngleRMS += trackCkvAngle*trackCkvAngle;

		runRingRadius += trackRadius;
		runTrackRadiusRMS += trackRadius*trackRadius;

		hTrackCkvAngle->Fill(trackCkvAngle/mrad);
		hTrackRadius->Fill(trackRadius/mm);

		for(Int_t i=0; i<nLayers; i++)
			hLayerNumPe[i]->Fill(layerNumPe[i]);

		if( betaDataExist ) {
			trackBeta = findMaxLikelihood(NumPhotoelectrons,vRadii);
			for(Int_t i=0; i<nLayers; i++)
			{
				if (layerNumPe[i]>0) {
					Double_t layerBeta = findMaxLikelihood(layerNumPe[i],vLayerRadii[i],i);
					hLayerTrackBeta[i]->Fill(layerBeta);
				}
			}
			runBeta += trackBeta;
			runTrackBetaRMS += trackBeta*trackBeta;
			hTrackBeta->Fill(trackBeta);
		}
	}

	hNumPhotoelectrons->Fill(NumPhotoelectrons);

	if (verboseLevel>1) {
		ancout << " Number of photoelectrons " << NumPhotoelectrons << "\n"
		       << " Cherenkov angle per track "
			   << trackCkvAngle << " mrad\n"
			   << " Ring radius per track " << trackRadius << " mm\n"
			   << " Beta per track " << trackBeta << G4endl;
	}
}


