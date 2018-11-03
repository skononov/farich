#include "RichDetectorConstruction.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <libgen.h>

#include <G4Box.hh>
#include <G4ThreeVector.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4SDManager.hh>
#include <G4VisAttributes.hh>
#include <G4ProcessTable.hh>

#include <G4GeometryManager.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>

#include <G4UnitsTable.hh>
#include <G4Element.hh>
#include <G4ElementTable.hh>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <G4MaterialTable.hh>
#include <G4MaterialPropertyVector.hh>

#include "RichGlobals.hh"
#include "RichDetectorMessenger.hh"
#include "RichParameters.hh"
#include "RichSpectrum.hh"
#include "RichAerogelMaterial.hh"
#include "MLRadiatorDescription.hh"
#include "RichPMT.hh"
#include "RichCerenkov.hh"
#include "RichFieldSetup.hh"

using namespace std;

RichDetectorConstruction::RichDetectorConstruction(int verb) :
	G4VUserDetectorConstruction()
{
	worldPhys			=0;
	radiatorPhys		=0;
	pmtPhys				=0;
	radiatorLog			=0;
	fieldSetup			=0;
	pmtSD				=0;
	runManager			=RichRunManager::GetRunManager();
	RichOpPhCathSurfacePT=0;
	RichQuartzMPT       =0;
	RichLiFMPT          =0;
	RichNaFMPT          =0;
	meanPhotonMomentum  =kReferencePhotMom; //400nm
	verboseLevel        =verb;
	detectorMessenger	=new RichDetectorMessenger(this);
	mlrDesc 			=new MLRadiatorDescription(proximityDistance,betaOptimized);

	SetDefaults();

	DefineMaterials();
}

RichDetectorConstruction::~RichDetectorConstruction()
{
//	G4MaterialTable::const_iterator cur=G4Material::GetMaterialTable()->begin(),
//			last=G4Material::GetMaterialTable()->end();

//	for ( ; cur!=last; cur++ ) delete *cur;

	delete RichAirMPT;
	delete RichAerogelMat;
	delete RichOpPhCathSurfacePT;
	delete RichQuartzMPT;
	delete RichLiFMPT;
	delete RichNaFMPT;
	delete mlrDesc;

	delete detectorMessenger;
}

void RichDetectorConstruction::DefineMaterials()
{
	G4double a,z;  //a=mass of a mole; z=mean number of protons;
	G4double density;
	G4String name,symbol;
	G4double defValue;

	G4int numel,natoms;  //numel=Number of elements constituting a material.

	G4UnitDefinition::BuildUnitsTable();

	G4NistManager* nistManager = G4NistManager::Instance();

	G4cout << "\nDefine Materials" << G4endl;

    G4cout << *(G4Element::GetElementTable()) << G4endl;

//Hydrogen
	a=1.00795*g/mole;
	G4Element* elH = new G4Element(name="Hydrogen",symbol="H",z=1.,a);

//Lithium
	a=6.941*g/mole;
	G4Element* elLi = new G4Element(name="Lithium",symbol="Li",z=3.,a);

//Nitrogen
	a=14.0067*g/mole;
	G4Element* elN = new G4Element(name="Nitrogen",symbol="N", z=7., a);

//Oxygen
	a=15.9994*g/mole;
	G4Element* elO = new G4Element(name="Oxygen",symbol="O", z=8., a);

//Fluorine
	a=18.9984*g/mole;
	G4Element* elF = new G4Element(name="Fluorine",symbol="F", z=9., a);

//Sodium
	a=22.9898*g/mole;
	G4Element* elNa = new G4Element(name="Sodium",symbol="Na", z=11., a);

//Magnesium
	a=24.305*g/mole;
	G4Element* elMg = new G4Element(name="Magnesium",symbol="Mg", z=12., a);

//Silicon
	a=28.0855*g/mole;
	G4Element* elSi = new G4Element(name="Silicon",symbol="Si",z=14.,a);

//Calcium
	a=40.078*g/mole;
	G4Element* elCa = new G4Element(name="Calcium",symbol="Ca",z=20.,a);

//Default optical properties
	G4double SimplePhotMomVector[2] = { kPhotonMinEnergy, kPhotonMaxEnergy };
	G4double DefAbsorpLength[2]  = { 1.0e32*mm, 1.0e32*mm };
	G4double DefRefIndex[2]      = { 1.0, 1.0 };

//Air at 20 degree C and 1 atm for the ambiet air.
	G4cout << " Defining Air .." << G4endl;
	G4Material* Air	= nistManager->FindOrBuildMaterial("G4_AIR");

	RichAirMPT = new G4MaterialPropertiesTable();

	RichAirMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);
	RichAirMPT->AddProperty("RINDEX",SimplePhotMomVector,DefRefIndex,2);

	Air->SetMaterialPropertiesTable(RichAirMPT);

//Water
	G4cout << " Defining Water .." << G4endl;
	density=1.000*g/cm3;
//	G4Material* H2O = new G4Material(name="Water",density,numel=2);
//	H2O->AddElement(elH,natoms=2);
//	H2O->AddElement(elO,natoms=1);
	G4Material* H2O = nistManager->FindOrBuildMaterial("G4_WATER");

	G4cout << "  Reading H2O refraction index from " << kH2OrefIndDataFile << G4endl;
	RichSpectrum riH2Osp(kH2OrefIndDataFile);
	G4cout << "   " << riH2Osp.Size() << " entries read" << G4endl;
	riH2Osp.ReduceToRange(kPhotonMinWL/nm,kPhotonMaxWL/nm);
	if (verboseLevel>0) riH2Osp.Print();

	RichH2OMPT = new G4MaterialPropertiesTable;
	RichH2OMPT->AddProperty("RINDEX",riH2Osp.GetMPV());
	RichH2OMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	H2O->SetMaterialPropertiesTable(RichH2OMPT);
    DefaultMPTs[H2O] = RichH2OMPT;

//SiO2
	G4cout << " Defining Fused Quartz .." << G4endl;
	density=2.200*g/cm3;
//	G4Material* SiO2 = new G4Material(name="Quartz",density,numel=2);
//	SiO2->AddElement(elSi,natoms=1);
//	SiO2->AddElement(elO,natoms=2);
    G4Material* SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	G4MaterialPropertyVector* quartzRefIndData = new G4MaterialPropertyVector;
	for (size_t ibin=0; ibin<defPhotMomVector.size(); ibin++) {
		G4double PhotMom = defPhotMomVector[ibin];
		quartzRefIndData->AddElement(PhotMom,FusedSilicaRefIndex(PhotMom));
	}

	RichQuartzMPT = new G4MaterialPropertiesTable;
	RichQuartzMPT->AddProperty("RINDEX",quartzRefIndData);
	RichQuartzMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	SiO2->SetMaterialPropertiesTable(RichQuartzMPT);
    DefaultMPTs[SiO2] = RichQuartzMPT;

//Aerogels
	G4cout << " Defining Aerogels .." << G4endl;
	RichAerogelMat = new RichAerogelMaterial(aerogelNomScatLength, aerogelAbsLenDataFile,
		chromaticity);
/*	density=0.25*g/cm3; //exact density doesn't matter because it is ideal simulation

	char strname[20];

    //predefine kMaxNumberOfLayers of possible aerogels
	for (G4int i=0; i<kMaxNumberOfLayers; i++) {
		sprintf(strname,"Aerogel%d",i+1);
		G4Material* anAerogel = new G4Material(G4String(strname), density, 2);
		anAerogel->AddMaterial(SiO2, 97.0*perCent);
		anAerogel->AddMaterial(H2O, 3.0*perCent); //small fraction of water
//		layerParam->AddAerogel(anAerogel);
		aerogelMaterials.push_back(anAerogel);
	}
*/

//Lithium fluoride (LiF)
	G4cout << " Defining Lithium Fluoride .." << G4endl;
	density=2.632*g/cm3; //PDG 2004

//	G4Material* LiF = new G4Material(name="LiF",density,numel=2);
//	LiF->AddElement(elLi,natoms=1);
//	LiF->AddElement(elF ,natoms=1);
	G4Material* LiF = nistManager->FindOrBuildMaterial("G4_LITHIUM_FLUORIDE");

	G4cout << "  Reading LiF refraction index from " << kLiFrefIndDataFile << G4endl;
	RichSpectrum riLiFsp(kLiFrefIndDataFile);
	G4cout << "   " << riLiFsp.Size() << " entries read" << G4endl;
	riLiFsp.ReduceToRange(kPhotonMinWL/nm,kPhotonMaxWL/nm);
	if (verboseLevel>0) riLiFsp.Print();

	RichLiFMPT = new G4MaterialPropertiesTable;
	RichLiFMPT->AddProperty("RINDEX",riLiFsp.GetMPV());
	RichLiFMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	LiF->SetMaterialPropertiesTable(RichLiFMPT);
    DefaultMPTs[LiF] = RichLiFMPT;

//Sodium fluoride (NaF)
	G4cout << " Defining Sodium Fluoride .." << G4endl;
	density=2.558*g/cm3; //PDG 2004

	G4Material* NaF = new G4Material(name="NaF",density,numel=2);
	NaF->AddElement(elNa,natoms=1);
	NaF->AddElement(elF ,natoms=1);

	G4cout << "  Reading NaF refraction index from " << kNaFrefIndDataFile << G4endl;
	RichSpectrum riNaFsp(kNaFrefIndDataFile);
	G4cout << "   " << riNaFsp.Size() << " entries read" << G4endl;
	riNaFsp.ReduceToRange(kPhotonMinWL/nm,kPhotonMaxWL/nm);
	if (verboseLevel>0) riNaFsp.Print();

	RichNaFMPT = new G4MaterialPropertiesTable;
	RichNaFMPT->AddProperty("RINDEX",riNaFsp.GetMPV());
	RichNaFMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	NaF->SetMaterialPropertiesTable(RichNaFMPT);
    DefaultMPTs[NaF] = RichNaFMPT;

//Calcium fluoride (CaF2)
	G4cout << " Defining Calcium Fluoride .." << G4endl;

	G4Material* CaF2 = nistManager->FindOrBuildMaterial("G4_CALCIUM_FLUORIDE");

	G4cout << "  Reading CaF2 refraction index from " << kCaF2refIndDataFile << G4endl;
	RichSpectrum riCaF2sp(kCaF2refIndDataFile);
	G4cout << "   " << riCaF2sp.Size() << " entries read" << G4endl;
	riCaF2sp.ReduceToRange(kPhotonMinWL/nm,kPhotonMaxWL/nm);
	if (verboseLevel>0) riCaF2sp.Print();

	RichCaF2MPT = new G4MaterialPropertiesTable;
	RichCaF2MPT->AddProperty("RINDEX",riCaF2sp.GetMPV());
	RichCaF2MPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	CaF2->SetMaterialPropertiesTable(RichCaF2MPT);
    DefaultMPTs[CaF2] = RichCaF2MPT;

//Plexiglas (PMMA)
	G4cout << " Defining Plexiglas .." << G4endl;
	G4Material* PMMA = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");

	G4cout << "  Reading Plexiglas refraction index from " << kPMMArefIndDataFile << G4endl;
	RichSpectrum riPMMAsp(kPMMArefIndDataFile);
	G4cout << "   " << riPMMAsp.Size() << " entries read" << G4endl;
	riPMMAsp.ReduceToRange(kPhotonMinWL/nm,kPhotonMaxWL/nm);
	if (verboseLevel>0) riPMMAsp.Print();

	RichPMMAMPT = new G4MaterialPropertiesTable;
	RichPMMAMPT->AddProperty("RINDEX",riPMMAsp.GetMPV());
	RichPMMAMPT->AddProperty("ABSLENGTH",SimplePhotMomVector,DefAbsorpLength,2);

	PMMA->SetMaterialPropertiesTable(RichPMMAMPT);
    DefaultMPTs[PMMA] = RichPMMAMPT;

//PMT material (vacuum)
	G4cout << " Defining PMT material .." << G4endl;

	new G4Material(name="PMTMaterial",z=1.,a=1.01*g/mole,
		density=universe_mean_density,kStateGas,0.1*kelvin,1e-19*pascal);

	//
	// Now for the material properties of Surfaces
	//
	G4cout << " Defining optical surfaces .." << G4endl;

	RichOpPhotocathodeSurface = new G4OpticalSurface("RichOpPhCathSurface",
		glisur,polishedbackpainted,dielectric_dielectric); //photons do not come out of PC

}

G4VPhysicalVolume* RichDetectorConstruction::Construct()
{
	G4cout << "\n*************BEGIN DETECTOR CONSTRUCTION*************" << G4endl << G4endl;

	ConstructDetector();

	G4cout << "\n**************END DETECTOR CONSTRUCTION**************" << G4endl << G4endl;

	fieldSetup = new RichFieldSetup;

	return worldPhys;
}

G4VPhysicalVolume* RichDetectorConstruction::ConstructDetector()
{
	if (worldPhys) {
		// clean-up previous geometry
		G4GeometryManager::GetInstance()->OpenGeometry();

		G4PhysicalVolumeStore::GetInstance()->Clean();
		G4LogicalVolumeStore::GetInstance()->Clean();
		G4SolidStore::GetInstance()->Clean();
		G4LogicalSkinSurface::CleanSurfaceTable();
		G4LogicalBorderSurface::CleanSurfaceTable();
		radiatorPhys=0;
		pmtPhys=0;
		worldPhys=0;
	}

//Optical properties of photocathode.
	DefineOpPhCathProperties();

//Create radiator logical volume.
//Uses mean photon momentum defined in DefineOpPhCathProperties().
	ConstructRadiator();

// World volume size along Z
//	G4double worldHalfZ=(proximityDistance+radiatorThickness+kDetectorThickness)/2;

// Radiator size & position along Z
	G4double radHalfZ=radiatorThickness/2;
	radiatorZposition=-proximityDistance/2-radHalfZ;

// PMT layer size&position along Z
	G4double detHalfZ=kDetectorThickness/2;
	G4double detPosZ=proximityDistance/2+detHalfZ;

//World volume
	G4VSolid *worldBox = new G4Box("WorldBox",kWorldHalfX,kWorldHalfY,
															kWorldHalfZ);
	G4LogicalVolume *worldLog =
		new G4LogicalVolume(worldBox,G4Material::GetMaterial("G4_AIR"),"WorldLogical");
	worldPhys =
		new G4PVPlacement(0,G4ThreeVector(0,0,0),worldLog,
			"WorldPhysical",0,0,0);

//Radiator placement
	if (radiatorLog)
		radiatorPhys = new G4PVPlacement(0,G4ThreeVector(0,0,radiatorZposition),radiatorLog,
			"RadiatorPhysical",worldLog,0,0);

//PMT layer
	G4VSolid *PMTLayer =
		new G4Box("PMTLayer",kDetHalfX,kDetHalfY,detHalfZ);
	G4LogicalVolume *PMTLog =
		new G4LogicalVolume(PMTLayer,G4Material::GetMaterial("PMTMaterial"),"PMTLogical");
	pmtPhys =
		new G4PVPlacement(0,G4ThreeVector(0,0,detPosZ),PMTLog,"PMTPhysical",worldLog,0,0);

//Photocathode surface
	new G4LogicalSkinSurface(G4String("PhotocathodeSurface"),
		PMTLog,RichOpPhotocathodeSurface);

//Now for the sensitive Detector
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	if( !pmtSD ) {
		pmtSD = new RichPMT("/RichPMT");
		SDman->AddNewDetector(pmtSD);
	}
	pixelSpacing=pixelSize/sqrt(geomEfficiency);
	pmtSD->SetPixelSizeAndSpacing(pixelSize,pixelSpacing,roundPixel);

//	PMTLog->SetSensitiveDetector(pmtSD);

//Now for visual attributes...
	G4VisAttributes* worldVisAttr =
		new G4VisAttributes(true,G4Colour(0.8,0.8,0.8));
	worldVisAttr->SetLineStyle(G4VisAttributes::dashed);
	worldVisAttr->SetForceWireframe(true);
	worldLog->SetVisAttributes(worldVisAttr);

	if( radiatorLog ) {
		radiatorLog->SetVisAttributes(new G4VisAttributes);
		for (G4int iDaughter=0; iDaughter<radiatorLog->GetNoDaughters(); iDaughter++) {
			G4String matName=radiatorLog->GetDaughter(iDaughter)->
				GetLogicalVolume()->GetMaterial()->GetName();
			G4Colour color;
			G4VisAttributes* layerVisAttr;
			if( matName.contains("Aerogel") ) {
                G4double f=0.05/(mlrDesc->GetIndex(iDaughter)-1);
				color=G4Colour(f*0.6,f*0.7,f*1.0,0.3);
			} else if( matName=="NaF" ) {
				color=G4Colour(1.0,0.0,1.0,0.5);
			} else if( matName=="LiF" ) {
				color=G4Colour(0.8,0.0,0.8,0.5);
			} else if( matName=="Water" ) {
				color=G4Colour(0.2,0.6,1.0,0.7);
			} else if( matName=="Quartz" || matName=="G4_PLEXIGLASS" ) {
				color=G4Colour(0.8,0.8,0.8,0.7);
			} else {
				color=G4Colour(1.0,1.0,0.8,0.5);
			}
			layerVisAttr=new G4VisAttributes(true,color);
			layerVisAttr->SetForceSolid(true);
			layerVisAttr->SetForceAuxEdgeVisible(true);
            radiatorLog->GetDaughter(iDaughter)->GetLogicalVolume()->
				SetVisAttributes(layerVisAttr);
		}
	}

	G4VisAttributes* PMTVisAttr =
		new G4VisAttributes(true,G4Colour(0.52,0.39,0.29));
	PMTVisAttr->SetForceSolid(true);
	PMTVisAttr->SetForceAuxEdgeVisible(true);
	PMTLog->SetVisAttributes(PMTVisAttr);

	modified = false;

//Print detector parameters
	Print();

	return worldPhys;
}

G4LogicalVolume* RichDetectorConstruction::ConstructRadiator()
{
	G4cout << " Constructing radiator in " << ModeName(mlrMode) << " mode .." << G4endl;

	RichAerogelMat->SetRefScatLength(aerogelNomScatLength);
	RichAerogelMat->SetAbsDataFile(aerogelAbsLenDataFile);

	RichAerogelMat->Initialize(verboseLevel);

	SwitchChromaticity();
	RichAerogelMat->SetChromaticity(chromaticity);
	if (verboseLevel) {
		G4cout << "  Chromaticity model set to ";
		if( chromaticity==kNoDispersion )
			G4cout << "'no disprersion'";
        else if( chromaticity==kLHCbModel )
			G4cout << "'LHCb model'";
        else
			G4cout << "'Quartz model'";
        G4cout << G4endl;
	}

	if (verboseLevel && mlrMode&fixedModes)
		G4cout << "  Set optimization beta to " << betaOptimized << G4endl;

	if (mlrMode==fixedNT)
        proximityDistance=overallDimension-totalThickness;

	if (verboseLevel)
		G4cout << "  Set proximity distance to " << proximityDistance/mm << " mm" << G4endl;

	mlrDesc->SetProximityDistance(proximityDistance);
	G4int N=0;

	if (mlrMode==manual) {
		N = mlrDesc->GetNlayers();
	} else if (mlrMode==fixedThickness) {
		if (N > kMaxNumberOfLayers) {
			G4cerr << "Number of layers too large. Maximum "
				   << kMaxNumberOfLayers << " allowed." << G4endl;
			return 0;
		}
		if (overallDimension < proximityDistance) {
			G4cerr << "Illegal parameters combination. "
				      "Dimension less than proximity distance." << G4endl;
			return 0;
		}
		mlrDesc->SetBeta(betaOptimized);
		G4double riPC=AerogelRefIndex(layer1NomRefIndex,meanPhotonMomentum);
		N = mlrDesc->MakeGabarit(overallDimension,riPC,normThickness1);
	} else if (mlrMode==fixedNumberOfLayers) {
		if (nLayers > kMaxNumberOfLayers) {
			G4cerr << "Number of layers too large. Maximum "
				   << kMaxNumberOfLayers << " allowed." << G4endl;
			return 0;
		}
		mlrDesc->SetBeta(betaOptimized);
		G4double riPC=AerogelRefIndex(layer1NomRefIndex,meanPhotonMomentum);
		N = mlrDesc->MakeLayers(nLayers,riPC,layer1Thickness);
		if (N>0 && N < nLayers) {
			G4cerr << "Failed to make " << nLayers << " layers. Proceeding with "
				   << N << " layers." << G4endl;
		}
	} else if (mlrMode==fixedNT) {
		if (nLayers > kMaxNumberOfLayers) {
			G4cerr << "Number of layers too large. Maximum "
				   << kMaxNumberOfLayers << " allowed." << G4endl;
			return 0;
		}
		if (overallDimension < totalThickness) {
			G4cerr << "Illegal parameters combination. "
				      "Dimension less than radiator thickness." << G4endl;
			return 0;
		}
		mlrDesc->SetBeta(betaOptimized);
		G4double riPC=AerogelRefIndex(layer1NomRefIndex,meanPhotonMomentum);
		N = mlrDesc->MakeFixed(nLayers,overallDimension,riPC);
		if (N < nLayers) {
			G4cerr << "Failed to make " << nLayers << " layers in radiator of"
				   << " total thickness " << totalThickness/mm << " mm" << G4endl;
			mlrDesc->clear();
			return 0;
		}
	} else { //mlrMode==multiRing
		if (nLayers > kMaxNumberOfLayers) {
			G4cerr << "Number of layers too large. Maximum "
				   << kMaxNumberOfLayers << " allowed." << G4endl;
			return 0;
		}
		G4double ri1PC=AerogelRefIndex(layer1NomRefIndex,meanPhotonMomentum);
		G4double ri2PC=AerogelRefIndex(series2NomRefIndex,meanPhotonMomentum);
		N = mlrDesc->MakeMultiRing(nLayers,overallDimension,ri1PC,ri2PC);
	}

	if( thicknessDistortion!=0.0 ) {
		G4cout << "  Distort thickness of layers by " << thicknessDistortion/perCent << "% "
			   << (evenDistortion?"evenly":"oddly") << G4endl;
		mlrDesc->DistortThickness(thicknessDistortion,evenDistortion);
	}

	if( densityDistortion!=0.0 ) {
		G4cout << "  Distort n-1 of layers by " << densityDistortion/perCent << "% "
			   << (evenDistortion?"evenly":"oddly") << G4endl;
		mlrDesc->DistortIndex(densityDistortion,evenDistortion);
	}

	if( uniformityDistortion!=0.0 ) {
		G4cout << "  Distort uniformity of layers by " << uniformityDistortion/perCent << "% "
			   << (evenDistortion?"evenly":"oddly") << G4endl;
		N = mlrDesc->DistortUniformity(uniformityDistortion,evenDistortion);
		if (N > kMaxNumberOfLayers) {
			G4cerr << "Uniformity distortion: Number of layers too large. Maximum "
				   << kMaxNumberOfLayers << " allowed." << G4endl;
			return 0;
		}
	}

	nLayers = N;

	radiatorThickness = mlrDesc->GetTotalThickness(); //overall radiator thickness

	if (!nLayers || radiatorThickness==0.0) {
		G4cerr << "No valid radiator defined" << G4endl;
		return 0;
	}

	if (proximityDistance/2+radiatorThickness > kWorldHalfZ ||
		proximityDistance/2+kDetectorThickness > kWorldHalfZ) {
		G4cerr << "Illegal geometry defined."
			   << " Either PMT or radiator come out of world volume." << G4endl;
		return 0;
	}

	if (verboseLevel) {
		G4cout << "  Number of layers in radiator: " << nLayers << G4endl;
		G4cout.setf(ios::fixed|ios::left);
		for (G4int i=0; i<nLayers; i++) {
			G4cout << "   Layer " << setw(2) << i+1 << ": "
				   << setw(10) << mlrDesc->GetMaterialName(i);
		if (G4String(mlrDesc->GetMaterialName(i)).contains("Aerogel"))
				G4cout << " n_nom=" << setprecision(4) << setw(6) << mlrDesc->GetIndex(i) << ",";
			G4cout << " t=" << setprecision(2) << mlrDesc->GetThickness(i)/mm << " mm"
				   << G4endl;
		}
		G4cout.unsetf(ios::fixed|ios::left);
		G4cout << "  Total thickness of radiator: " << radiatorThickness/mm << " mm" << G4endl;
		G4cout.precision(0);
		G4cout.width(0);
	}

// Radiator mother volume - "Holder"
	G4VSolid *radiatorBox = new G4Box("RadiatorHolder",
		kRadHalfX,kRadHalfY,radiatorThickness/2);

	G4Material* Air = G4Material::GetMaterial("G4_AIR");

	radiatorLog =
		new G4LogicalVolume(radiatorBox,Air,"RadiatorHolder");

	G4double z = radiatorThickness/2; //z-position of the downstream face of radiator
	char strname[20];
	RichUsedMaterials.resize(nLayers);

// Assign properties to aerogels, construct solids and logical volumes of layers
	for (G4int i=0; i<nLayers; i++) {
		G4double ri=0, halfT, halfPitch;
		G4Material* aMaterial=0;
		RichUsedMaterials[i]=Air;

		halfT = halfPitch = mlrDesc->GetThickness(i)/2;
		G4String matName(mlrDesc->GetMaterialName(i));
		if (matName=="Aerogel") {
			//set mean refractive index in MLRadiatorDescription for manual mode
			if (mlrMode==manual)
				ri = AerogelRefIndex(mlrDesc->GetIndex(i),meanPhotonMomentum);
			else
				ri = mlrDesc->GetIndex(i);
			if (ri==1.0) { //skip this layer, to be filled by air
				z -= 2*halfPitch;
				continue;
			}
			aMaterial = RichAerogelMat->GetAerogelWithIndex(ri,meanPhotonMomentum);
		} else if (matName=="Air") {
            //just skip this layer, to be filled by air
			z -= 2*halfPitch;
			continue;
		} else {
			halfT -= 0.005*mm; //make a small air gap between layers of materials
			aMaterial = G4Material::GetMaterial(matName);
			if( !aMaterial ) {
				G4cerr << "Undefined material " << matName << " met. Make a gap." << G4endl;
				z -= 2*halfPitch;
				continue;
			}
			if (!aMaterial->GetMaterialPropertiesTable()) {
				G4cerr << "Material " << matName << " does not have MPT. Make a gap." << G4endl;
				z -= 2*halfPitch;
				continue;
			}
			G4MaterialPropertyVector* riV=aMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
			if (!riV) {
				G4cerr << "Material " << matName << " does not have RINDEX property. Make a gap." << G4endl;
				z -= 2*halfPitch;
				continue;
			}
            //fill in refractive index for this layer in MLRadiatorDescription
			mlrDesc->SetIndex(i,riV->GetProperty(meanPhotonMomentum));
		}

		RichUsedMaterials[i]=aMaterial;

		sprintf(strname,"Layer%d",i+1);

		G4VSolid *layerBox = new G4Box(strname,kRadHalfX,kRadHalfY,halfT);
		G4LogicalVolume *layerLog = new G4LogicalVolume(layerBox,aMaterial,strname);

		z -= halfPitch;

		new G4PVPlacement(0,G4ThreeVector(0,0,z),layerLog,strname,radiatorLog,0,i);

		z -= halfPitch;
	}

	//set NULL property table for every unused aerogel materials
//	for (G4int i=aerogelIndex; i<kMaxNumberOfLayers; i++)
//		aerogelMaterials[i]->SetMaterialPropertiesTable(0);

	return radiatorLog;
}

void RichDetectorConstruction::AddLayer(G4String layerSpec,G4double thickness)
{
	Changed(manual);

	char* endptr;
	const char* nstr=layerSpec.data();

	G4double number = strtod(nstr,&endptr);

//	G4cout<<"spec='"<<layerSpec<<"' number="<<number<<" endptr="<<(void*)endptr<<G4endl;

	if( *endptr==0 ) { // number
		if ( number<1.0 || number==-HUGE_VAL || number==HUGE_VAL ) {
			G4cerr << "Invalid refractive index value " << number << " ignored" << G4endl;
			return;
		}
		if ( number==1.0 ) {
			mlrDesc->AddAlayer("Air",thickness);
			G4cout << "Gap added with size of " << thickness/mm << " mm" <<G4endl;
		} else {
			mlrDesc->AddAlayer(number,thickness);
			G4cout << "Aerogel layer added with nominal index of refraction of " << number
				   << " and thickness of " << thickness/mm << " mm" <<G4endl;
		}
	} else { // string
		if ( layerSpec.compareTo("lif",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("LiF",thickness);
			G4cout << "LiF layer added of thickness "
				   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("naf",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("NaF",thickness);
			G4cout << "NaF layer added of thickness "
				   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("caf2",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("G4_CALCIUM_FLUORIDE",thickness);
			G4cout << "CaF2 layer added of thickness "
				   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("pmma",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("G4_PLEXIGLASS",thickness);
			G4cout << "G4_PLEXIGLASS layer added of thickness "
				   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("quartz",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("Quartz",thickness);
			G4cout << "Quartz layer added of thickness "
				   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("water",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("Water",thickness);
			G4cout << "Water layer added of thickness "
                                   << thickness/mm << " mm" <<G4endl;
		} else if ( layerSpec.compareTo("air",G4String::ignoreCase)==0 ||
			        layerSpec.compareTo("none",G4String::ignoreCase)==0 ) {
			mlrDesc->AddAlayer("Air",thickness);
			G4cout << "Gap added with size of " << thickness/mm << " mm" <<G4endl;
		} else {
			G4cerr << "Unknown material name " << layerSpec << " ignored" << G4endl;
			return;
		}
	}
}

void RichDetectorConstruction::DefineOpPhCathProperties()
{
	G4cout << " Defining optical photocathode properties .." << G4endl;

	if (RichOpPhCathSurfacePT) {
		RichOpPhCathSurfacePT->RemoveProperty("RINDEX");
		RichOpPhCathSurfacePT->RemoveProperty("REFLECTIVITY");
		RichOpPhCathSurfacePT->RemoveProperty("EFFICIENCY");
	}
	else
		RichOpPhCathSurfacePT = new G4MaterialPropertiesTable;

    if (verboseLevel)
		G4cout << "  Reading quantum efficiency from " << qeDataFile << G4endl;
	RichSpectrum pmtQEsp(qeDataFile);

	if (pmtQEsp.IsEmpty()) {
    	if (verboseLevel)
			G4cout << "   Set quantum efficiency to 1.0" << G4endl;
		pmtQEsp.Set(kPhotonMinWL/nm,kPhotonMaxWL/nm,1.0);
	} else {
		if (verboseLevel) {
			G4cout << "   " << pmtQEsp.Size() << " entries read" << G4endl;
			if (verboseLevel>1) pmtQEsp.Print();
		}
	}

    if (verboseLevel)
		G4cout << "  Set detector acceptance to " << detectionEfficieny/perCent << "%" << G4endl;

	G4MaterialPropertyVector *pcEffVector = new G4MaterialPropertyVector;
	G4MaterialPropertyVector *pcReflVector = new G4MaterialPropertyVector;
	G4MaterialPropertyVector *pcRefIndVector =
		new G4MaterialPropertyVector(*RichQuartzMPT->GetProperty("RINDEX"));

	RichSpectrumIter next(pmtQEsp);
	pair<double,double> *pp;

	if (verboseLevel>1)
		G4cout << "  Photodetection efficiency:\n  eV   %" << G4endl;

	G4double PhotMom=0;
	G4double minPhotWL, maxPhotWL;
	pmtQEsp.GetRange(minPhotWL,maxPhotWL);

	//calculate mean photon momentum for this photocathode
	G4double norm=0, prevPhotMom;
	meanPhotonMomentum=0;

	while( (pp=next()) ) {
		prevPhotMom=PhotMom;
		PhotMom = kPhotMomWaveConv/pp->first*eV;
		G4double qri = pcRefIndVector->GetProperty(PhotMom);
		G4double eff = detectionEfficieny*(pp->second/100.)*(qri+1)*(qri+1)/4/qri;
		//G4cout<<PhotMom/eV<<"\t"<<eff<<G4endl;
		pcEffVector->AddElement(PhotMom,eff);
		if( prevPhotMom!=0 ) {
			meanPhotonMomentum+=eff*PhotMom*(prevPhotMom-PhotMom);
			norm+=eff*(prevPhotMom-PhotMom);
		}
		if (verboseLevel>1)
			G4cout<<setw(4)<<setprecision(2)<<PhotMom/eV<<"  "<<eff/perCent<<G4endl;
	}

	meanPhotonMomentum/=norm;

	if (verboseLevel) {
		G4cout << "  Mean detected photon momentum: " << setw(5) << setprecision(4)
			   << meanPhotonMomentum/eV << " eV" << " ("
			   << kPhotMomWaveConv/(meanPhotonMomentum/eV) << "nm wavelength)" << G4endl;

		G4cout.precision(0);
		G4cout.width(0);
	}

	G4double minPhotMom=(kPhotMomWaveConv/maxPhotWL-0.001)*eV,
	         maxPhotMom=(kPhotMomWaveConv/minPhotWL+0.001)*eV;
	pcEffVector->AddElement(kPhotonMinEnergy,0.0);
	pcEffVector->AddElement(minPhotMom,0.0); //force zero efficiency at the ends
	pcEffVector->AddElement(maxPhotMom,0.0); //of spectrum range
	pcEffVector->AddElement(kPhotonMaxEnergy,0.0);
	pcReflVector->AddElement(kPhotonMinEnergy,0.0);
	pcReflVector->AddElement(3.0*eV,0.0);
	pcReflVector->AddElement(kPhotonMaxEnergy,0.0);

	RichOpPhCathSurfacePT->AddProperty("EFFICIENCY",pcEffVector);
	RichOpPhCathSurfacePT->AddProperty("REFLECTIVITY",pcReflVector);
	RichOpPhCathSurfacePT->AddProperty("RINDEX",pcRefIndVector);

	RichOpPhotocathodeSurface->
		SetMaterialPropertiesTable(RichOpPhCathSurfacePT);
}

void RichDetectorConstruction::SwitchChromaticity()
{
    MMPT_t::const_iterator first = DefaultMPTs.begin(), last = DefaultMPTs.end();
	for( ; first!=last; first++) {
	    G4Material* mat = first->first;
	    G4MaterialPropertiesTable* def_mpt = first->second;
		if (mat->GetMaterialPropertiesTable()!=def_mpt)
			delete mat->GetMaterialPropertiesTable();
		if (chromaticity==kNoDispersion) {
			G4MaterialPropertiesTable* tmp_mpt = new G4MaterialPropertiesTable;
			G4MaterialPropertyVector* riVector = new G4MaterialPropertyVector;
			G4double ri = def_mpt->GetProperty("RINDEX")->GetProperty(meanPhotonMomentum);
			riVector->AddElement(kPhotonMinEnergy,ri);
			riVector->AddElement(kPhotonMaxEnergy,ri);
			tmp_mpt->AddProperty("RINDEX",riVector);
			tmp_mpt->AddProperty("ABSLENGTH",
				new G4MaterialPropertyVector(*def_mpt->GetProperty("ABSLENGTH")));
			mat->SetMaterialPropertiesTable(tmp_mpt);
		} else {
			mat->SetMaterialPropertiesTable(def_mpt);
		}
	}
	if (chromaticity==kNoDispersion)
		AerogelRefIndex = AerogelRefIndexConstant;
	else if(chromaticity==kLHCbModel)
		AerogelRefIndex = AerogelRefIndexLHCb;
	else
		AerogelRefIndex = AerogelRefIndexQuartz;
}

void RichDetectorConstruction::Update()
{
	if (!modified) return;

	G4cout << "\n*************BEGIN DETECTOR UPDATE*************" << G4endl << G4endl;

	// clean-up previous geometry
	G4GeometryManager::GetInstance()->OpenGeometry();

	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();
	G4LogicalSkinSurface::CleanSurfaceTable();
	G4LogicalBorderSurface::CleanSurfaceTable();

	// define new one
	worldPhys = ConstructDetector();

	if( !worldPhys ) {
		G4cerr << " Construction of the geometry failed! " << G4endl;
		return;
	}

	runManager->DefineWorldVolume(worldPhys,true);

//	runManager->GeometryHasBeenModified();
	runManager->PhysicsHasBeenModified();

	G4VProcess* cp=G4ProcessTable::GetProcessTable()->FindProcess("Cerenkov","e-");
	if (cp) ((RichCerenkov*)cp)->RebuildPhysicsTable();
//	if (cp) ((RichCerenkov*)cp)->ResetPhysicsTable();

	G4cout << "\n**************END DETECTOR UPDATE**************" << G4endl << G4endl;
}

void RichDetectorConstruction::SetDefaults()
{
	mlrDesc->clear();

	mlrMode=fixedNumberOfLayers;

	layer1NomRefIndex=defAerogelNominalRefIndex;
	proximityDistance=defProximityDistance;
	overallDimension=defOverallDimension;
	betaOptimized=1.0;
	nLayers=0;
	layer1Thickness=0;
	normThickness1=0;

	//no distortion by default
	thicknessDistortion =0.0;
	densityDistortion   =0.0;
	uniformityDistortion=0.0;
	evenDistortion = true;

	chromaticity = kLHCbModel;

	aerogelNomScatLength=defAerogelScatLength;
	aerogelBoundaryLoss=defAerogelSurfTrans;
	aerogelAbsLenDataFile=defAerogelAbsLenDataFile;

	detectionEfficieny=defDetectionEfficiency;
	geomEfficiency=defGeomEfficiency;
	pixelSize=pixelSpacing=0.0;
	roundPixel=false;
	qeDataFile=defPMTQEdataFile;

	Changed();
}

void RichDetectorConstruction::Print() const
{
	char* fnAbs = basename((char*)aerogelAbsLenDataFile.c_str());
	char* fnQE  = basename((char*)qeDataFile.c_str());

	G4cout << "-----------------------------------------------\n"
		   << "|             Detector parameters             |\n"
		   << "-----------------------------------------------\n"
		   << "Aerogel:\n"
		   << "  Scattering length: " << aerogelNomScatLength/mm << " mm\n"
		   << "  Absorption data file: " << fnAbs << "\n"
		   << "Radiator:\n"
           << "  Constructed in " << ModeName(mlrMode) << " mode\n"
		   << "  Number of layers: " << nLayers << G4endl;
	G4cout.setf(ios::fixed|ios::left);
	for (G4int i=0; i<nLayers; i++) {
		G4cout << "  Layer " << setw(2) << i+1 << ": ";
		G4Material* aMaterial=RichUsedMaterials[i];
		if( !aMaterial ) {
			G4cerr << " Can not find material for the layer" << G4endl;
			break;
		}
		if (aMaterial->GetName()=="G4_AIR")
			G4cout << setw(10) << "Gap";
		else {
			G4cout << setw(10) << aMaterial->GetName()
				   << " rho=" << setprecision(4) << aMaterial->GetDensity()/g*cm3 << " g/cm^3,";
			G4double riRef=1.0;
			G4MaterialPropertiesTable* mpt=aMaterial->GetMaterialPropertiesTable();
			if (mpt) {
				G4MaterialPropertyVector* riV=mpt->GetProperty("RINDEX");
                if (riV)
					riRef=riV->GetProperty(kReferencePhotMom);
			}
			G4cout << " n(400nm)=" << setprecision(4) << riRef << ",";
		}
		G4cout << " t=" << setprecision(2) << mlrDesc->GetThickness(i)/mm << "mm" << G4endl;
	}
	G4cout.unsetf(ios::fixed|ios::left);
	G4cout.width(0);

	G4cout << "  Total thickness: " << setprecision(4) << radiatorThickness/mm << " mm\n"
	       << "Beta to optimize radiator: " << setprecision(6) << betaOptimized << "\n"
		   << "Proximity distance: " << setprecision(4) << proximityDistance/mm << " mm\n"
		   << "Photodetector:\n"
		   << "  QE data file: " << fnQE << "\n"
		   << "  detection efficiency: " << detectionEfficieny/perCent << "%\n"
		   << "  geometrical efficiency: " << geomEfficiency/perCent << "%\n";
	if( pixelSize>0.0 )
		G4cout << "  pixel size: " << pixelSize/mm << "mm\n";
	else
		G4cout << "  absolute position resolution\n";
	G4cout << "  mean photon wavelength: " << kPhotMomWaveConv/(meanPhotonMomentum/eV) << " nm"
		   << G4endl;
	G4cout.precision(0);

	if (modified)
		G4cout << "Need update!" << G4endl;
	G4cout << "-----------------------------------------------" << G4endl;
}

G4double RichDetectorConstruction::GetNpeRoughly(G4double beta,G4double incAngle)
{
	const G4double PhotMomBin = defPhotMomVector[1]-defPhotMomVector[0];
	const G4double CherenkovSpectrumConst=fine_structure_const/hbarc;

	G4MaterialPropertyVector* qeData =
		RichOpPhCathSurfacePT->GetProperty("EFFICIENCY");

	vector<G4MaterialPropertyVector*> riV(nLayers);
	vector<G4MaterialPropertyVector*> absLengthV(nLayers);
	vector<G4MaterialPropertyVector*> scatLengthV(nLayers);

	for (G4int i=0; i<nLayers; i++) {
		riV[i]=RichUsedMaterials[i]->
			GetMaterialPropertiesTable()->GetProperty("RINDEX");
		absLengthV[i]=RichUsedMaterials[i]->
			GetMaterialPropertiesTable()->GetProperty("ABSLENGTH");
		scatLengthV[i]=RichUsedMaterials[i]->
			GetMaterialPropertiesTable()->GetProperty("RAYLEIGH");
	}

	vector<G4double> ri(nLayers), att(nLayers);

	G4double integral=0;

	for (size_t ibin=0; ibin<defPhotMomVector.size()-1; ibin++) {
		G4double PhotMom = defPhotMomVector[ibin];
        G4double eff = geomEfficiency*qeData->GetProperty(PhotMom);

		G4double spectralDensity=0;
		for (G4int i=0; i<nLayers; i++) {
			ri[i]=riV[i]->GetProperty(PhotMom);

			G4double cosT = 1/beta/ri[i];

			if (cosT>1) continue; //below threshold

			if (ri[i]*sin(acos(cosT)-incAngle)>=1.0) continue; //total internal reflection

			G4double labs=absLengthV[i]?absLengthV[i]->GetProperty(PhotMom):1E+32*mm;
			G4double lsc=scatLengthV[i]?scatLengthV[i]->GetProperty(PhotMom):1E+32*mm;

			att[i]=1/labs+1/lsc;

			G4double I=CherenkovSpectrumConst*eff*(1-cosT*cosT);

			G4double path=mlrDesc->GetThickness(i)/cos(incAngle);

			if (att[i]*path<1e-15) I*=path;
			else I*=(1-exp(-att[i]*path))/att[i];

			//photon angle to layer's normal (if there is incAngle least angle is taken)
			G4double initAngleCos=cos(acos(cosT)-incAngle);

			for (G4int j=i-1; j>=0; j--) {
				if (ri[j]==1.0) continue; //no attenuation
				G4double photonPath=mlrDesc->GetThickness(j) /
					sqrt(1-ri[i]*ri[i]/(ri[j]*ri[j])*(1-initAngleCos*initAngleCos));
				I*=exp(-att[j]*photonPath);
			}

			spectralDensity += I;
		}
		integral += spectralDensity * PhotMomBin;
	}

	return integral;
}

