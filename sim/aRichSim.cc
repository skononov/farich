#include <iostream>
#include <unistd.h>
#include "RichRunAction.hh"
#include "RichEventAction.hh"
#include "RichDetectorConstruction.hh"
#include "RichPrimaryGeneratorAction.hh"
#include "RichTrackingAction.hh"
#include "RichSteppingAction.hh"
#include "RichPhysicsList.hh"
//#include "RichAnalysisManagerLH.hh"
#include "RichAnalysisManagerS.hh"
#include "RichRunManager.hh"

#include <Randomize.hh>
#include <G4VPhysicalVolume.hh>
#include <G4UImanager.hh>
#include <G4UIterminal.hh>

#ifdef G4UI_USE_XM
# include <G4UIXm.hh>
#endif
#ifdef G4VIS_USE
# include "RichVisManager.hh"
#endif
#ifdef G4UI_USE_TCSH
# include <G4UItcsh.hh>
#endif

#include <G4ios.hh>
#include <cstdlib>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;

static const char *progname="aRichSim";
static const char *optstring="hvb:";
static int detverblev=0;
static bool batch=false;
static const char *macro=0;

static void Usage(int status)
{
	cout<<"Usage:\n"
		<<progname<<" [-v] [-b macro | macro]\n"
		<<"Geant4 FARICH simulation\n"
		<<" macro - macro file to execute before the interactive session starts\n"
		<<"Options:\n"
		<<"  -v               Some detailed printout in RichDetectorConstruction\n"
		<<"  -v -v            More detailed printout in RichDetectorConstruction\n"
		<<"  -b macro         Batch mode, execute macro and exit.\n"
		<<"  -h               Print out this help\n"
		<<endl;
	exit(status);
}

int main(int argc, char **argv)
{
// Check program options
	int opt;

	while( (opt=getopt(argc,argv,optstring))>0 ) {
		switch( opt ) {
			case '?': Usage(1); break;
			case 'h': Usage(0); break;
			case 'v': detverblev++; break;
			case 'b': batch=true; macro=optarg; break;
			default : cerr<<"Missing option "<<opt<<endl; Usage(1);
		}
	}

	if( optind<argc ) {
		if( batch ) {
			cerr<<"Macro file is already given for the batch mode"<<endl;
			Usage(1);
		}
		macro=argv[optind];
		if( optind+1<argc )
			cerr<<"Only the first argument accepted as a macro name. The rest is ignored."<<endl;
	}

	if( macro && access(macro,R_OK)!=0 ) {
		cerr<<"Can not read from the file "<<macro<<endl;
		return 1;
	}


// Choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

// Seed the random number generator manually
    time_t curtime = time(NULL);
	G4long myseed = curtime;
	CLHEP::HepRandom::setTheSeed(myseed);

	cout<<ctime(&curtime)<<endl;

// Run manager

	RichRunManager *runManager = new RichRunManager;

	RichDetectorConstruction *RichDet = new RichDetectorConstruction(detverblev);
	runManager->SetUserInitialization(RichDet);

	RichPhysicsList *RichPL = new RichPhysicsList;
	runManager->SetUserInitialization(RichPL);

// UserAction classes - optional

#ifdef G4VIS_USE
	RichVisManager *visManager = new RichVisManager;
	visManager->SetVerboseLevel(0);
	visManager->Initialize();
#endif

	RichVAnalysisManager *analysis = RichAnalysisManagerS::getInstance();

	runManager->SetUserAction(new RichRunAction(analysis));

	runManager->SetUserAction(new RichPrimaryGeneratorAction(RichDet));

	runManager->SetUserAction(new RichEventAction(analysis));

	runManager->SetUserAction(new RichTrackingAction);

	runManager->SetUserAction(new RichSteppingAction);

//	RichStackingAction *stackAction = new RichStackingAction;
//	runManager->SetUserAction(stackAction);

	G4UImanager *UI = G4UImanager::GetUIpointer();

	G4UIsession *session = 0;

	//Initialize G4 kernel
	runManager->Initialize();

    if( macro )
		UI->ExecuteMacroFile(macro);

	if ( !batch ) { // Interactive mode
#ifdef G4UI_USE_TCSH
		session = new G4UIterminal(new G4UItcsh);
#else
		session = new G4UIterminal();
#endif
		session->SessionStart();
		delete session;
	}

#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;

	delete analysis;

	return 0;
}
