#include <fstream>
#include <cctype>
#include <cmath>

#include <G4Timer.hh>
#include <G4ios.hh>
#include <G4Run.hh>
#include <G4HCofThisEvent.hh>

//ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TTree.h"

#include "RichVAnalysisManager.hh"
#include "RichGlobals.hh"
#include "RichHit.hh"
#include "RichAnalysisMessenger.hh"
#include "RichFieldSetup.hh"
#include "RichDetectorConstruction.hh"
#include "RichPrimaryGeneratorAction.hh"
#include "RichPrimaryData.hh"
#include "RichRootEvent.hh"
#include "RichParticlePropagation.hh"
#include "MLRadiatorDescription.hh"

#include "G4SystemOfUnits.hh"

#define SUM2_to_RMS(sum2,mean,n) (sum2<n*mean*mean?0:sqrt(sum2/n-mean*mean))

RichVAnalysisManager *RichVAnalysisManager::instance = 0;
RichDetectorConstruction *RichVAnalysisManager::detector = 0;
RichPrimaryGeneratorAction *RichVAnalysisManager::primary = 0;

TROOT root("root","Rich analysis");

RichVAnalysisManager::RichVAnalysisManager() :
	ancout(G4cout.rdbuf()),
	currentRun(0),
	outputFile(0),
	outputFileName("rich.root"),
	macroFileName("run.C"),
	betaDataFileName(),
	interactive(true),
	verboseLevel(1),
	saveHitTree(true),
	saveEventTree(false),
	weighPixels(false),
	betaDataExist(false)
{
	iTimer = new G4Timer;

	int argc=2;
	char *argv[]={"rint","-l"};
	application = new TRint("rint",&argc,argv);

	analysisMessenger = new RichAnalysisMessenger(this);

    rootEvent = new RichRootEvent;

	ancout.precision(5);
	ancout.width(7);
}

RichVAnalysisManager::~RichVAnalysisManager()
{
	delete application;
	delete analysisMessenger;
	if(outputFile) delete outputFile;
	delete iTimer;
	delete rootEvent;
	RichRootEvent::Reset();
	instance = 0;
}

void RichVAnalysisManager::book(const G4Run* run)
{
	currentRun = (G4Run*)run;

	if (!primary)
		primary = (RichPrimaryGeneratorAction*)
			G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

	if (!detector)
		detector = (RichDetectorConstruction*)
			G4RunManager::GetRunManager()->GetUserDetectorConstruction();

	beta = primary->GetBeta();
//	G4double bg = beta/sqrt(1-beta*beta);

	//estimates for position and direction of the primary
	primVertexPosition = primary->GetPosition();
 	primDirection = primary->GetMomentumDirection();
	momentum = primary->GetTotalMomentum();
	charge = primary->GetCharge();

	//get magnetic field value
	magField = detector->GetFieldSetup()->GetFieldValue();

	ancout << " ----------------- Input data -----------------\n"
		   << " Primary particle: " << primary->GetParticleName() << "\n"
		   << "           charge: " << int(charge) << " e\n"
		   << "         momentum: " << momentum/GeV << " GeV\n"
		   << "             beta: " << beta << "\n"
		   << "       beta*gamma: " << 1/sqrt(1/beta/beta-1) << "\n"
		   << "  vertex position: " << primVertexPosition/mm << " mm\n"
		   << "        direction: " << primDirection <<", theta="
		   << primDirection.theta()/degree << " deg, phi="
		   << primDirection.phi()/degree << " deg\n"
		   << "   magnetic field: " << magField/tesla << " T\n";

	//estimate of mean emission point of optical photons
//	meanEmissionPoint = position +
//		direction * (detector->GetRadiatorZposition()-position.z())/direction.z();

	RichParticlePropagation pp(magField,primVertexPosition,momentum*primDirection,charge);

	G4ThreeVector position, direction, dummy;
	G4bool res = pp.Propagate(detector->GetRadiatorZposition(),position,direction);
	if( res ) {
		meanEmissionPoint = position;
		meanPrimDirection = direction;
	} else {
		position = meanEmissionPoint = G4ThreeVector(0.,0.,detector->GetRadiatorZposition());
		direction = meanPrimDirection = G4ThreeVector(0.,0.,1.);
	}

	eventVertexShift = G4ThreeVector(0,0,0);

//	inPrimaryPosition = position +
//		direction * (detector->GetRadiatorInSurfZ()-position.z())/direction.z();
//	outPrimaryPosition = position +
//		direction * (detector->GetRadiatorOutSurfZ()-position.z())/direction.z();

	//estimate of primary crossing positions for input and output radiator surfaces
	pp.Propagate(detector->GetRadiatorInSurfZ(),inPrimaryPosition,dummy);
	pp.Propagate(detector->GetRadiatorOutSurfZ(),outPrimaryPosition,dummy);

	ancout << "  mean vertex pos: " << meanEmissionPoint/mm << " mm\n"
		   << "   mean direction: " << meanPrimDirection <<", theta="
		   << meanPrimDirection.theta()/degree << " deg, phi="
		   << meanPrimDirection.phi()/degree << " deg" << G4endl;

	incAngle = direction.theta(); //particle direction angle from normal
	G4double tgIncAngle=tan(incAngle);

	nLayers = detector->GetNumberOfLayers();
	mlrMode = detector->GetMLRmode();
	nRings = mlrMode==RichDetectorConstruction::multiRing?nLayers:1;

	proximityDistance = detector->GetProximityDistance();
	MLRadiatorDescription* mlrDesc=detector->GetMLRadiatorDescription();

	G4double z=detector->GetRadiatorOutSurfZ();

	//find maximum ring radius and determine Cherenkov cone angle for layers
	G4double maxRingRadius=0;
	G4double maxSpotRadius=0;
	G4double maxAngle=0;
	for (G4int i=0; i<nLayers; i++) {
		layerIndices[i]=1.02*mlrDesc->GetIndex(i); //take into account dispersion
		layerThicknesses[i]=mlrDesc->GetThickness(i);
		layerAngles[i]=0;

		z-=0.5*layerThicknesses[i];
		pp.Propagate(z,layerVertexPositions[i],dummy);
		z-=0.5*layerThicknesses[i];

		if (beta < 1/layerIndices[i]) continue;

		layerAngles[i]=acos(1/beta/layerIndices[i]);
		if (maxAngle<layerAngles[i]) maxAngle=layerAngles[i];

		if (layerAngles[i]+incAngle>halfpi ||
			layerIndices[i]*sin(layerAngles[i]+incAngle)>1) {
			maxRingRadius=HUGE_VAL;
			break;
		}

		G4double tgi=tan(layerAngles[i]+incAngle);

		G4double radius=(tgi-tgIncAngle)*layerThicknesses[i];
		for(G4int j=i-1; j>=0; j--)
			radius+=(tanrefangle(tgi,layerIndices[i],layerIndices[j])-tgIncAngle)*layerThicknesses[j];

		if (radius>maxSpotRadius) maxSpotRadius=radius;

		radius+=(tanrefangle(tgi,layerIndices[i],1.0)-tgIncAngle)*proximityDistance;

		if (radius>maxRingRadius) maxRingRadius=radius;
	}

    //Take values of radiator index for the mean wavelength
	for (G4int i=0; i<nLayers; i++)
		layerIndices[i]=mlrDesc->GetIndex(i);

	imageCenter = position +
		direction*(detector->GetDetectorWindowZ()-position.z())/direction.z();

	G4double maxXimageHalfSize=fabs(copysign(kDetHalfX,imageCenter.x())-imageCenter.x());
	G4double maxYimageHalfSize=fabs(copysign(kDetHalfY,imageCenter.y())-imageCenter.y());
	G4double imageHalfSize=maxXimageHalfSize>maxYimageHalfSize?maxXimageHalfSize:maxYimageHalfSize;

	if (isinf(maxRingRadius)) {
		maxRingRadius = imageHalfSize;
		ancout << " Warning: Total internal reflection is possible! The ring will be stripped." << G4endl;
	} else {
		if( fabs(copysign(maxRingRadius,imageCenter.x())+imageCenter.x()) > kDetHalfX ||
			fabs(copysign(maxRingRadius,imageCenter.y())+imageCenter.y()) > kDetHalfY )
			ancout << " Warning: Part of the ring can be stripped by photodetector boundaries!" << G4endl;
		else {
			if( maxRingRadius*2+10*mm < imageHalfSize )
				imageHalfSize = maxRingRadius*2+10*mm;
		}
	}

	maxAngle=maxAngle*1.2+100*mrad;

	expectedNpe = detector->GetNpeRoughly(beta,incAngle);
    G4int npeRange = int( expectedNpe + 5*sqrt(expectedNpe) + 2 );

	ancout << " Roughly estimated number of photoelectrons: " << expectedNpe << "\n"
		   << " Max spot radius at radiator: " << maxSpotRadius/mm << " mm\n"
		   << " Max ring radius: " << maxRingRadius/mm << " mm\n"
		   << " Average ring center: " << imageCenter/mm << " mm\n"
		   << " Weigh pixels by number of photoelectrons: " << (weighPixels?"true":"false") << "\n"
		   << " Store hit tree: " << (saveHitTree?"true":"false") << "\n"
		   << " Store event tree: " << (saveEventTree?"true":"false") << "\n"
		   << "-----------------------------------------------" << G4endl;

	//Load beta vs radius precalculated values if any
	LoadBetaData();

// Book histograms
	outputFile = new TFile(outputFileName,"RECREATE","FARICH simulation results",5);

	if( saveHitTree ) {
		tHits = new TTree("th","Hits information");
		tHits->Branch("hit",&hitinfo,hitinfofmt);
        tHits->SetAutoSave(10000000);
		tHits->SetMarkerStyle(1);
	} else
		tHits = 0;

	if( saveEventTree ) {
		tEvent = new TTree("te","Event information");
        tEvent->SetAutoSave(20000000);
		tEvent->Branch("event",&rootEvent,16000,2);
		tEvent->SetMarkerStyle(1);
	} else
		tEvent = 0;

	//Store primary data to the tree. Note: no primary vertex position is stored as it varies from event to event
	RichPrimaryData primdata(primary->GetParticle(),
		primary->GetTotalMomentum()/GeV, direction.theta()/degree, direction.phi()/degree);

	primdata.Write();

	if( saveEventTree ) {
		rootEvent->particle.name = (const char*)primary->GetParticleName();
		rootEvent->particle.m = primary->GetMass()/GeV;
	}

	hNumCkvPhot =
		new TH1F("hnphot","Number of Cherenkov photons per track;N",500,0.,1000.);
	hNumCkvPhot->SetBit(TH1::kCanRebin);
	hNumPhotAtPMT =
		new TH1F("hnphotd","Number of photons at PMT;N",500,0.,500.);
	hNumPhotAtPMT->SetBit(TH1::kCanRebin);
	hNumPhotoelectrons =
		new TH1F("hnpe","Number of photoelectrons per track;N_{pe}",npeRange+1,-.5,npeRange+0.5);
	hPeWL =
		new TH1F("hpewl","Wavelength distribution of photons to hit;wavelength, nm",
			83,170.,1000.);
	hNumPcIncidenceVsAngle =
		new TH1F("hNumPcInc","Distribution of photons at PC on angle of incidence;"
		"angle of incidence, deg",200,0.,90.);
	hPcReflectionVsAngle =
		new TH1F("hPcRefl","Reflecion coefficient of PC vs angle of incidence;"
		"angle of incidence, deg",200,0.,90.);

	hXYhits = new TH2F("hxy","XY hit map at PMT;X, mm;Y, mm",
		1000,(imageCenter.x()-imageHalfSize)/mm,(imageCenter.x()+imageHalfSize)/mm,
		1000,(imageCenter.y()-imageHalfSize)/mm,(imageCenter.y()+imageHalfSize)/mm);
	hXYhits->SetMarkerStyle(1);

	hCkvAngle = new TH1F("hang","Cherenkov angle per photon; #theta_{c}, mrad",
		1000,0,maxAngle/mrad);

	hAlphaAngle = new TH1F("halp","Azimutal angle per photon; #alpha, degrees",
		200,-180,180);

	pCkvAngVsAlpha = new TProfile("pangalp","Cherenkov angle vs azimuthal angle;#alpha, degrees;#theta_{c}, mrad;",
		60,-180,180,"s");

	h1peRadius = new TH1F("hrad","Cerenkov ring radius per photon;radius, mm",
		int(imageHalfSize/(0.5*mm)),0,imageHalfSize/mm);

	h1peBeta = new TH1F("hbet","#beta distribution per photon;#beta",1000,beta-.2,beta+.2);

	hTrackCkvAngle = new TH1F("hang_t","Cherenkov angle per track;#theta_{c}, mrad",1500,
		0,maxAngle/mrad);

	hTrackRadius = new TH1F("hrad_t","Cerenkov ring radius per track;radius, mm",int(imageHalfSize/(0.1*mm)),
		0,imageHalfSize/mm);

	hTrackBeta = new TH1F("hbet_t","Measured #beta distribution;#beta",1500,beta-.2,beta+.2);

	hNpePerPixel = new TH1F("hnpe_px","Number of photoelectrons per pixel;N_{pe}",10,0,10);

	char str[10], title[100];
	for (G4int i=0; i<nLayers; i++) {
		sprintf(str,"hnpe%d",i+1);
		hLayerNumPe[i] = (TH1F*)hNumPhotoelectrons->Clone(str);
		sprintf(title,"Number of photoelectrons for layer %d",i+1);
		hLayerNumPe[i]->SetTitle(title);
		hLayerNumPe[i]->SetLineColor(i+1);

		sprintf(str,"hang%d",i+1);
		hLayerCkvAngle[i] = (TH1F*)hCkvAngle->Clone(str);
		sprintf(title,"Cerenkov angle for layer %d",i+1);
		hLayerCkvAngle[i]->SetTitle(title);
		hLayerCkvAngle[i]->SetLineColor(i+1);

		sprintf(str,"hrad%d",i+1);
		hLayer1peRadius[i] = (TH1F*)h1peRadius->Clone(str);
		sprintf(title,"Cerenkov ring radius for layer %d",i+1);
		hLayer1peRadius[i]->SetTitle(title);
		hLayer1peRadius[i]->SetLineColor(i+1);

		sprintf(str,"hbet%d",i+1);
		hLayer1peBeta[i] = (TH1F*)h1peBeta->Clone(str);
		sprintf(title,"#beta distribution for layer %d per photon",i+1);
		hLayer1peBeta[i]->SetTitle(title);
		hLayer1peBeta[i]->SetLineColor(i+1);

		sprintf(str,"hbet%d_t",i+1);
		hLayerTrackBeta[i] = (TH1F*)hTrackBeta->Clone(str);
		sprintf(title,"#beta distribution for layer %d per track",i+1);
		hLayerTrackBeta[i]->SetTitle(title);
		hLayerTrackBeta[i]->SetLineColor(i+1);
	}
}

void RichVAnalysisManager::finish()
{
	hPcReflectionVsAngle->Divide(hNumPcIncidenceVsAngle);

	ancout << "Storing histograms in " << outputFileName << "..." << G4endl;
	outputFile->Write();
	if (interactive) {
		if (gSystem->AccessPathName(macroFileName)==0) {
			ancout << "Executing root macro " << macroFileName << "..." << G4endl;
			gROOT->Macro(macroFileName);
		}
		application->Run(kTRUE); //run CINT at the end of a run
	}

	outputFile->Close();
	delete outputFile;
	outputFile=0;
}

void RichVAnalysisManager::BeginOfEventAnalysis(const G4Event *event)
{
	iTimer->Start();
	NumCkvPhot = 0;
	NumPhotoelectrons = 0;
	NumPhotAtPMT = 0;

	//primary position for this event
	primVertexPosition = primary->GetPosition();

    //shift of primary vertex in XY-plane w.r.t. X=Y=0 line for the event
	eventVertexShift = primVertexPosition - G4ThreeVector(0.,0.,primVertexPosition.z());

	//initialize primary direction
	primDirection = primary->GetMomentumDirection();

	if( saveEventTree ) {
		rootEvent->Clear();
		rootEvent->particle.x0 = meanEmissionPoint.x()+eventVertexShift.x();
		rootEvent->particle.y0 = meanEmissionPoint.y()+eventVertexShift.y();
		rootEvent->particle.z0 = meanEmissionPoint.z()+eventVertexShift.z();
		rootEvent->particle.vx = meanPrimDirection.x();
		rootEvent->particle.vy = meanPrimDirection.y();
		rootEvent->particle.vz = meanPrimDirection.z();
		rootEvent->particle.p = primary->GetTotalMomentum()/GeV;
	}

	//initialize in and out positions
//	inPrimaryPosition = primVertexPosition +
//		primDirection * (detector->GetRadiatorInSurfZ()-primVertexPosition.z())/primDirection.z();
//	outPrimaryPosition = primVertexPosition +
//		primDirection * (detector->GetRadiatorOutSurfZ()-primVertexPosition.z())/primDirection.z();

	if (verboseLevel>1) {
		ancout << "Event " << event->GetEventID() << " analysis:" << G4endl;
	}
	else if (verboseLevel==1 && (event->GetEventID()%10)==0) {
		G4cerr << "Event " << event->GetEventID() << "           \r"
			   << std::flush;
	}
}

void RichVAnalysisManager::find_thetac(const G4ThreeVector& hitpos,G4int iLayer,G4double& thetac,G4double& alpha)
{
	G4double theta=0;
	G4ThreeVector w;

	if (nRings>1) {
		G4double D2=proximityDistance*proximityDistance;
		w = hitpos-layerVertexPositions[iLayer]-eventVertexShift;

		G4double x=w.perp();
		G4double n2=layerIndices[iLayer]*layerIndices[iLayer];

		G4double Tn=tan(w.theta());
		G4double T=0;
		G4double xr;

		// iteration loop
		while (fabs(1-T/Tn)>1e-6) {
			T = Tn;
			xr = x - 0.5*T*layerThicknesses[iLayer];
			for (G4int j=iLayer-1; j>=0; j--)
				xr -= tanrefangle(T,layerIndices[iLayer],layerIndices[j])*layerThicknesses[j];
			Tn = xr/sqrt(D2*n2+xr*xr*(n2-1));
		}
		theta=atan(Tn);
	} else {
		w = hitpos-meanEmissionPoint-eventVertexShift;

		G4double a = w.perp()/detector->GetRadiatorThickness()*2;
		G4double b = w.z()/detector->GetRadiatorThickness()*2, c2 = (b-1)*(b-1);
		G4double n2 = layerIndices[0]*layerIndices[0];

		G4double Tn = tan(w.theta()); //zero approximation for T evalulation
		G4double T = 0;

		// iteration loop
		while (fabs(1-T/Tn)>1e-6) {
			T = Tn;
			Tn = 1 / sqrt(n2*(c2/(a-T)/(a-T)+1)-1);
		}
		theta=atan(Tn);
	}

    //Cherenkov photon direction in the lab frame
	G4ThreeVector photDir(1,0,0);
	photDir.setTheta(theta);
	photDir.setPhi(w.phi());

    photDir.rotateZ(-meanPrimDirection.phi());
    photDir.rotateY(-meanPrimDirection.theta());
	thetac = photDir.theta(); //calculated Cherenkov angle
	alpha =  photDir.phi();   //calculated azimuthal angle of Cherenkov photon
}

