#ifndef RichTrajectory_hh
# define RichTrajectory_hh 1

# include <G4Trajectory.hh>
# include <G4Allocator.hh>

class G4Polyline;
class G4Track;
class G4ParticleDefinition;


class RichTrajectory : public G4Trajectory {
public:
	RichTrajectory();
	RichTrajectory(const G4Track * aTrack);
	RichTrajectory(RichTrajectory &);

	virtual ~ RichTrajectory();

	virtual void DrawTrajectory(G4int i_mode = 0) const;

	inline void *operator  new(size_t);
	inline void operator  delete(void *);

	void SetDrawTrajectory(G4bool b) {	drawit = b;	}

private:
	G4bool drawit;
	G4ParticleDefinition *particleDefinition;
	G4double wavelength; //optical photon wavelength, nm
};

extern G4Allocator < RichTrajectory > RichTrajectoryAllocator;

inline void *RichTrajectory::operator  new(size_t)
{
	void *aTrajectory;
	aTrajectory = (void *) RichTrajectoryAllocator.MallocSingle();
	return aTrajectory;
}

inline void RichTrajectory::operator  delete(void *aTrajectory)
{
	RichTrajectoryAllocator.FreeSingle((RichTrajectory *) aTrajectory);
}

#endif
