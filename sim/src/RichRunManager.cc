#include "RichRunManager.hh"

RichRunManager* RichRunManager::fRunManager = 0;

RichRunManager* RichRunManager::GetRunManager()
{ return fRunManager; }

RichRunManager::RichRunManager() {
	if(fRunManager)
	{ G4Exception("RichRunManager constructed twice."); }
	fRunManager = this;
}


