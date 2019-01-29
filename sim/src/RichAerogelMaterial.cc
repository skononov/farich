#include "RichAerogelMaterial.hh"

#include <cmath>
#include <iomanip>
#include <libgen.h>

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "RichGlobals.hh"
#include "RichSpectrum.hh"

RichAerogelMaterial::RichAerogelMaterial(G4double Lsc,const G4String& AbsorpDataFile,
	G4int chrom) :
	refScatLength (Lsc),
	absorpDataFile (AbsorpDataFile),
	chromaticity (chrom),
	absorpData (0),
	rayleighData (0),
	initialized (false)
{
	Initialize();
}

G4int RichAerogelMaterial::Initialize(G4int verboseLevel)
{
// Load aerogel absorption length data
    if (verboseLevel)
		G4cout << "  Reading aerogel absorption length from "  << absorpDataFile << G4endl;

	RichSpectrum Lab_wl(absorpDataFile);

	absorpData = new G4MaterialPropertyVector;

	if (Lab_wl.IsEmpty()) {
		if (verboseLevel)
			G4cout << "   No absorption length data for aerogel." << G4endl;
		absorpData->InsertValues(kPhotonMinEnergy,1e32*CLHEP::mm);
		absorpData->InsertValues(kPhotonMaxEnergy,1e32*CLHEP::mm);
	} else {
		if (verboseLevel) {
			G4cout << "   " << Lab_wl.Size() << " entries read" << G4endl;
			if (verboseLevel>1) Lab_wl.Print();
		}
		Lab_wl.ReduceToRange(kPhotonMinWL/CLHEP::nm,kPhotonMaxWL/CLHEP::nm);

		RichSpectrumIter next(Lab_wl);
		std::pair<double,double> *pp;
		while( (pp=next()) ) {
			G4double PhotMom=kPhotMomWaveConv/pp->first*CLHEP::eV; //assuming WL in nm
			absorpData->InsertValues(PhotMom,pp->second*CLHEP::cm);
		}
	}

    if (verboseLevel)
		G4cout << "  Aerogel absorption length at " << kReferenceWavelength/CLHEP::nm
			   << " nm: " << absorpData->Value(kReferencePhotMom)/CLHEP::cm << " cm\n";

// Define Rayleigh scattering length vector
	rayleighData = new G4MaterialPropertyVector;
    if (verboseLevel)
		G4cout << "  Aerogel Rayleigh scattering length at " << kReferenceWavelength/CLHEP::nm
			   << " nm: " << refScatLength/CLHEP::cm << " cm\n";

	for (size_t ibin=0; ibin<defPhotMomVector.size(); ibin++) {
		G4double PhotMom = defPhotMomVector[ibin];
		rayleighData->InsertValues(PhotMom,refScatLength*pow(kReferencePhotMom/PhotMom,4));
	}

	initialized = true;

	return 0;
}

RichAerogelMaterial::~RichAerogelMaterial()
{
	if (absorpData) delete absorpData;
	if (rayleighData) delete rayleighData;
}

G4Material* RichAerogelMaterial::GetAerogelWithIndex(G4double ri,G4double pm)
{
	if (!initialized) Initialize();

	G4Material *anAerogel=0;

	G4MaterialPropertyVector* riData = new G4MaterialPropertyVector;

	G4double riRef = AerogelRefIndex(ri,pm,true);

	if (aerogels.find(riRef) != aerogels.end())
		anAerogel=aerogels[riRef];

	if (chromaticity==kNoDispersion) {
		riData->InsertValues(kPhotonMinEnergy,ri);
		riData->InsertValues(kPhotonMaxEnergy,ri);
	} else {
		for (size_t ibin=0; ibin<defPhotMomVector.size(); ibin++)
			riData->InsertValues(defPhotMomVector[ibin],AerogelRefIndex(riRef,defPhotMomVector[ibin],false));
	}

	if (!anAerogel) {
		G4double density = (riRef*riRef-1)/kAlpha400*CLHEP::g/CLHEP::cm3;

		char name[30];
		sprintf(name,"Aerogel_%6.4f",riRef);

        anAerogel = new G4Material(name, density, 2, kStateSolid);
//		anAerogel->AddMaterial(G4Material::GetMaterial("Quartz"), 97.0*perCent);
//		anAerogel->AddMaterial(G4Material::GetMaterial("Water"), 3.0*perCent); //small fraction of water
		anAerogel->AddMaterial(G4Material::GetMaterial("G4_SILICON_DIOXIDE"), 97.0*CLHEP::perCent);
		anAerogel->AddMaterial(G4Material::GetMaterial("G4_WATER"), 3.0*CLHEP::perCent); //small fraction of water

		aerogels[riRef] = anAerogel;
	}

	G4MaterialPropertyVector* absorpDataCopy = new G4MaterialPropertyVector(*absorpData);

	G4MaterialPropertyVector* rayleighDataCopy = new G4MaterialPropertyVector(*rayleighData);

	G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable;

	mpt->AddProperty("RINDEX", riData);
	mpt->AddProperty("ABSLENGTH", absorpDataCopy);
	mpt->AddProperty("RAYLEIGH", rayleighDataCopy);

	if (anAerogel->GetMaterialPropertiesTable())
		delete anAerogel->GetMaterialPropertiesTable();

	anAerogel->SetMaterialPropertiesTable(mpt);

	return anAerogel;
}

