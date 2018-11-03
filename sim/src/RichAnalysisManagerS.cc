#include <fstream>
#include <cctype>
#include <cmath>

#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4HCofThisEvent.hh"
#include "Randomize.hh"

//ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"

#include "RichGlobals.hh"
#include "RichHit.hh"
#include "RichAnalysisMessenger.hh"
#include "RichAnalysisManagerS.hh"
#include "RichDetectorConstruction.hh"
#include "RichPrimaryGeneratorAction.hh"
#include "RichRootEvent.hh"

#define SUM2_to_RMS(sum2,mean,n) (sum2<n*mean*mean?0:sqrt(sum2/n-mean*mean))

RichAnalysisManagerS::RichAnalysisManagerS() : RichVAnalysisManager()
{
	fBeta      = 0;
	fSigmaBeta = 0;
	fRing      = 0;
}

RichAnalysisManagerS::~RichAnalysisManagerS()
{}

void RichAnalysisManagerS::book(const G4Run* run)
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
	run1peBetaRMS=0;
	runTrackBetaRMS=0;

	run1peAngle=0;
	run1peRadius=0;

	runNumEvWithHit=0;
	runNumHits=0;
}

void RichAnalysisManagerS::finish()
{
// Some results
	runNumPe /= currentRun->GetNumberOfEvent();

	runCkvAngle      /= runNumEvWithHit;
	runRingRadius    /= runNumEvWithHit;
	run1peAngle      /= runNumHits;
	run1peRadius     /= runNumHits;
	run1peAngleRMS    = SUM2_to_RMS(run1peAngleRMS,run1peAngle,runNumHits);
	runTrackAngleRMS  = SUM2_to_RMS(runTrackAngleRMS,runCkvAngle,runNumEvWithHit);
	run1peRadiusRMS   = SUM2_to_RMS(run1peRadiusRMS,run1peRadius,runNumHits);
	runTrackRadiusRMS = SUM2_to_RMS(runTrackRadiusRMS,runRingRadius,runNumEvWithHit);
	runBeta         /= runNumEvWithHit;
	run1peBetaRMS   /= runNumHits;
	runTrackBetaRMS /= runNumEvWithHit;

	Int_t nwLayers=0;
	for (G4int i=0; i<nLayers; i++)
		if (hLayer1peRadius[i]->GetEntries()>0) nwLayers++;

	ancout << "------------------ Run " << currentRun->GetRunID() << " summary --------------------------" << "\n"
		   << "Number of events with at list one hit in this run " << runNumEvWithHit << "\n"
	 	   << " Number of working layers: " << nwLayers << "\n"
	 	   << " Mean number of photoelectrons: " << runNumPe << "\n"
		   << " Cherenkov angle, mrad: " << runCkvAngle/mrad << " measured, (1st layer "
		   << layerAngles[0]/mrad << " est)\n"
		   << " Cherenkov ring radius, mm: " << runRingRadius/mm << "\n"
		   << " Angular resolution, mrad: " << run1peAngleRMS/mrad << " per photon, "
		   << runTrackAngleRMS/mrad << " per track\n"
		   << " Radius resolution, mm: " << run1peRadiusRMS/mm << " per photon, "
		   << runTrackRadiusRMS/mm << " per track\n"
		   << " Velocity: " << runBeta << " (" << beta << "), diff=" << fabs(runBeta-beta) << "\n"
		   << " Velocity resolution: " << run1peBetaRMS << " per photon, "
		   << runTrackBetaRMS << " per track\n"
		   << "------------------------------------------------------------" << G4endl;

	RichVAnalysisManager::finish();
}

using namespace std;

#include "RichAnalysisModuleS.C"

void RichAnalysisManagerS::LoadBetaData()
{
	betaDataExist=false;

	if (fBeta)      delete fBeta;
	if (fSigmaBeta) delete fSigmaBeta;
	if (fRing)       delete fRing;
	fBeta      = 0;
	fSigmaBeta = 0;
	fRing       = 0;

	if (betaDataFileName.empty()) return;

	Int_t n=loadcal(betaDataFileName);

	if (!n) {
		G4cerr << "Calibration data incorrect and will not be used" << G4endl;
		return;
	}
	if (Nr!=nRings) {
		G4cerr << "Current number of rings does not match the number read from file" << G4endl;
		return;
	}
	if (incAngle==0 && Nphi>1 || incAngle>0 && Nphi==1) {
		G4cerr << "Particle angle does not match that read from file" << G4endl;
		return;
	}

	fBeta      = new TF1("fbeta",radius_to_beta,0,100,2);
	fBeta->SetParameters(0,0);
	fSigmaBeta = new TF1("fsigbeta",sigma_beta,0,1.1,2);
	fSigmaBeta->SetParameters(0,0);
	fRing       = new TF1("fring",ring_number,0,nRings,1);
	fRing->SetParameter(0,beta);

	betaDataExist=true;
}


void RichAnalysisManagerS::EndOfEventAnalysis(const RichHitsCollection *RHC)
{
	iTimer->Stop();

	if (verboseLevel>1) {
		ancout << *iTimer << "\n"
		       << " Number of Cherenkov photons " << NumCkvPhot << "\n"
		       << " Number of photons at PMT " << NumPhotAtPMT << G4endl;
	}

	if( hNumCkvPhot->GetXaxis()->GetXmax()<NumCkvPhot )
		hNumCkvPhot->RebinAxis(NumCkvPhot,hNumCkvPhot->GetXaxis());
	hNumCkvPhot->Fill(NumCkvPhot);
	if( hNumPhotAtPMT->GetXaxis()->GetXmax()<NumPhotAtPMT )
		hNumPhotAtPMT->RebinAxis(NumPhotAtPMT,hNumPhotAtPMT->GetXaxis());
	hNumPhotAtPMT->Fill(NumPhotAtPMT);

	if (RHC==0) {
        if( saveEventTree )
			tEvent->Fill();
		return;
	}

/*  Obsolete. Primary particle direction and position is now correctly calculated in RichVAnalysisManager
	G4ThreeVector diff=outPrimaryPosition-inPrimaryPosition;
//	if( diff.z()>0. ) {
		//recalculate primary particle parameters
		primDirection=diff.unit(); //estimate for mean primary direction in radiator
		meanEmissionPoint=0.5*(outPrimaryPosition+inPrimaryPosition);
		imageCenter=inPrimaryPosition +
			primDirection*(detector->GetDetectorWindowZ()-inPrimaryPosition.z())/primDirection.z();
		imageCenter = meanEmissionPoint +
			primMeanDirection*(detector->GetDetectorWindowZ()-meanEmissionPoint.z())/primMeanDirection.z();
		if (verboseLevel>1) {
			ancout << " Primary position\n  at input radiator surface: " << inPrimaryPosition/mm << "\n"
				   << "  at output radiator surface: " << outPrimaryPosition/mm << "\n"
				   << "  mean primary direction: " << primDirection << "\n"
				   << "  at photodetector plane: " << imageCenter/mm << G4endl;
		}
//	} else {
//		G4cerr << "\nInsane primary position at radiator boundaries! "
//			   << " Primary position at input radiator surface: " << inPrimaryPosition/mm << "\n"
//			   << "                  at output radiator surface: " << outPrimaryPosition/mm << G4endl;
//	}
*/

	G4double trackCkvAngle=0, trackRadius=0, trackBeta=0, trackSumOfWeights=0;
	G4double weight=1.0;

	memset(layerBeta,0,nLayers*sizeof(G4double));
	memset(layerNumPe,0,nLayers*sizeof(G4int));

	for (G4int ih=0; ih<RHC->entries(); ih++ )
	{
		RichHit* aHit = (RichHit*)RHC->GetHit(ih);

#ifndef BACKGROUND_ON
		if (!aHit->IsHit()) continue;
#endif

		if( weighPixels ) weight = aHit->GetNhits();

		NumPhotoelectrons += int(weight);

		G4ThreeVector hitPos = aHit->GetPosition();

		G4double x=hitPos.x()/mm, y=hitPos.y()/mm;

		hXYhits->Fill(x,y,weight);

		G4double wl_nm = aHit->GetWL()/nanometer;
		hPeWL->Fill(wl_nm,weight);

		hNpePerPixel->Fill(aHit->GetNhits());

		G4int iLayer = aHit->GetLayerNo();

		//Cherenkov angle fine estimation. Account for refraction.
		G4double thetac, alpha;
		find_thetac(hitPos,iLayer,thetac,alpha);
		if( TMath::IsNaN(thetac) ) thetac=0;

		//Ring radius calculation
		G4double radius = (hitPos-imageCenter-eventVertexShift).mag();

        //Get Cherenkov angle in mrad and radius in mm
		G4double thetac_mr = thetac/mrad, radius_mm=radius/mm;

		trackCkvAngle += thetac*weight;
		trackRadius += radius*weight;

		run1peAngle += thetac*weight;
		run1peRadius += radius*weight;

		run1peAngleRMS += thetac*thetac*weight;
		run1peRadiusRMS += radius*radius*weight;

		hCkvAngle->Fill(thetac_mr,weight);
		hAlphaAngle->Fill(alpha/degree,weight);
		h1peRadius->Fill(radius_mm,weight);
		pCkvAngVsAlpha->Fill(alpha/degree,thetac_mr,weight);

		hLayerCkvAngle[iLayer]->Fill(thetac_mr);
		hLayer1peRadius[iLayer]->Fill(radius_mm);
		layerNumPe[iLayer]++;

        //beta determination
		G4double b = 0.0;
		G4double sb = 1.0;

		G4int iRing=0;

		if (betaDataExist) {
			//determine a ring the photon belongs to
			iRing = (int)fRing->Eval(radius_mm);
			G4int iphi = int(acos(x/radius_mm)/dPhi);
			if (iphi>=Nphi) iphi=Nphi-1;

			fBeta->SetParameters(iRing,iphi);
			fSigmaBeta->SetParameters(iRing,iphi);

			//Beta calculation
			b = fBeta->Eval(radius_mm);
			sb = fSigmaBeta->Eval(b);
		}

		if (b>0) {
			G4double wb=1/(sb*sb)*weight;
			trackBeta += b*wb;
			layerBeta[iLayer] += b;
			trackSumOfWeights += wb;
			run1peBetaRMS += sb*weight;
			h1peBeta->Fill(b,weight);
            hLayer1peBeta[iLayer]->Fill(b);
		}

		//Copy hit information to the storable structures
		if( saveHitTree || saveEventTree ) {
			G4ThreeVector vertexPos(aHit->GetVertexPos());
			G4ThreeVector dir0(aHit->GetVertexDirection());
			G4ThreeVector dir0loc(aHit->GetVertexDirection());
            //transform dir0 to cf w.r.t. primary particle
			dir0loc.rotateZ(-primDirection.phi());
			dir0loc.rotateY(-primDirection.theta());
			G4ThreeVector dir(aHit->GetDirection());
			G4ThreeVector expos(aHit->GetExactPosition());

			if( saveHitTree ) {
				hitinfo.event=aHit->GetEvent();

				hitinfo.x=x;
				hitinfo.y=y;
				hitinfo.layer=iLayer+1;

				hitinfo.x0=vertexPos.x()/mm;
				hitinfo.y0=vertexPos.y()/mm;
				hitinfo.z0=vertexPos.z()/mm;

				hitinfo.theta0=dir0loc.theta()/degree;
				hitinfo.phi0=dir0loc.phi()/degree;

				hitinfo.theta=dir.theta()/degree;
				hitinfo.phi=dir.phi()/degree;

				hitinfo.wl=wl_nm;

				hitinfo.radius=radius_mm;
				hitinfo.angle=thetac_mr;
				hitinfo.alpha=alpha/degree;
				hitinfo.beta=b;

				hitinfo.xe=expos.x()/mm;
				hitinfo.ye=expos.y()/mm;

				hitinfo.npe=aHit->GetNhits();

				hitinfo.scattered=!aHit->IsHit();

				tHits->Fill();
			}

			if( saveEventTree ) {
				RichRootHit &rootHit = *rootEvent->AddHit();

				rootHit.x=x;
				rootHit.y=y;
        	    rootHit.z=detector->GetDetectorWindowZ()/mm;
				rootHit.layer=iLayer+1;

				rootHit.x0=vertexPos.x()/mm;
				rootHit.y0=vertexPos.y()/mm;
				rootHit.z0=vertexPos.z()/mm;

				rootHit.vx=dir0.x();
				rootHit.vy=dir0.y();
				rootHit.vz=dir0.z();
				rootHit.thetac0=dir0loc.theta();
				rootHit.phic0=dir0loc.phi();

				rootHit.wl=wl_nm;

				rootHit.radius=radius_mm;
				rootHit.thetac=thetac_mr/1e3;
				rootHit.phic=alpha;

				rootHit.xe=expos.x()/mm;
				rootHit.ye=expos.y()/mm;
				rootHit.ze=expos.z()/mm;

				rootHit.n=aHit->GetNhits();

				rootHit.scattered=!aHit->IsHit();
			}
		}

        //Print hit info
		if (verboseLevel>2) {
			ancout << "Hit#" << NumPhotoelectrons << "\n"
				   << " hitPos = " << hitPos/mm << "\n"
				   << " wavelength = " << wl_nm << "nm\n"
				   << " layer# = " << iLayer+1 << ", ring# = " << iRing+1 << "\n"
				   << " radius = " << radius_mm << "mm\n"
				   << " angle = " << thetac_mr << "mrad\n"
				   << " beta = " << b << "+-" << sb << G4endl;
		}
	}

	runNumPe += NumPhotoelectrons;
	hNumPhotoelectrons->Fill(NumPhotoelectrons);
    if( saveEventTree )
		tEvent->Fill();

	for(Int_t i=0; i<nLayers; i++)
	{
		hLayerNumPe[i]->Fill(layerNumPe[i]);
		if ( betaDataExist && layerNumPe[i]>0) {
			layerBeta[i] /= layerNumPe[i];
			hLayerTrackBeta[i]->Fill(layerBeta[i]);
		}
	}

	if (NumPhotoelectrons>0) {
		trackCkvAngle /= NumPhotoelectrons;
		trackRadius /= NumPhotoelectrons;

		runNumHits += NumPhotoelectrons;
		runNumEvWithHit++;

		runCkvAngle += trackCkvAngle;
		runTrackAngleRMS += trackCkvAngle*trackCkvAngle;

		runRingRadius += trackRadius;
		runTrackRadiusRMS += trackRadius*trackRadius;

		if (trackSumOfWeights>0) {
			G4double trackSigmaBeta = sqrt(1/trackSumOfWeights);
			trackBeta /= trackSumOfWeights;
			runBeta += trackBeta;
			runTrackBetaRMS += trackSigmaBeta;
			hTrackBeta->Fill(trackBeta);
		}

		hTrackCkvAngle->Fill(trackCkvAngle/mrad);
		hTrackRadius->Fill(trackRadius/mm);
	}

	if (verboseLevel>1) {
		ancout << " Number of photoelectrons " << NumPhotoelectrons << "\n"
		       << " Cherenkov angle per track "
			   << trackCkvAngle << " mrad\n"
			   << " Ring radius per track " << trackRadius << " mm\n"
			   << " Beta per track " << trackBeta << G4endl;
	}
}


