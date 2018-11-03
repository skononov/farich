#ifndef RichRunManager_h
#define RichRunManager_h 1

#include <G4RunManager.hh>

class RichRunManager : public G4RunManager
{
public: // with description
    static RichRunManager* GetRunManager();
    //  Static method which returns the singleton pointer of G4RunManager or
    // its derived class.

private:
    static RichRunManager* fRunManager;

public: // with description
	RichRunManager();

	~RichRunManager() {}

	inline void GeometryHasBeenModified()
	{
//		geometryNeedsToBeClosed = true;
		geometryInitialized = false;
//		kernel->GeometryHasBeenModified();
//		SetGeometryToBeOptimized(true);
	}
	inline void PhysicsHasBeenModified()
	{
//		physicsInitialized = false;
		kernel->PhysicsHasBeenModified();
	}
};

#endif
