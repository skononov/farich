#ifndef RichVisManager_h
# define RichVisManager_h 1

# ifdef G4VIS_USE

#  include <G4VisManager.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RichVisManager : public G4VisManager {
public:
	RichVisManager();
	virtual ~ RichVisManager();
private:
	void RegisterGraphicsSystems();
	void RegisterModelFactories();

};

# endif

#endif
