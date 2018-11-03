#ifndef RichHit_h
# define RichHit_h 1

# include <G4VHit.hh>
# include <G4Allocator.hh>
# include <G4ThreeVector.hh>
# include <G4THitsCollection.hh>

class RichVisManager;

class RichHit : public G4VHit {
public:
	RichHit() : G4VHit(),
		event(0),
		pixelId(0),
		pixelPosition(),
		exactPosition(),
		vertex(),
		direction(),
		wavelength(-1),
		nHit(0),
		layerNo(-1)
	{}
	virtual ~RichHit();

	G4int operator==(const RichHit &right) const
	{
		return pixelPosition==right.pixelPosition;
	}

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	void Draw();
	void Print();

private:
	G4int event; //event number where this hit appears
	G4int pixelId;
	G4ThreeVector pixelPosition; //hit position by center of pixel
	G4ThreeVector exactPosition; //hit exact position
	G4ThreeVector vertex; //origin of the photon
	G4ThreeVector vertexDirection, direction; //actual photon direction in vertex and at PMT entrance
	G4double wavelength; //photon wavelength
	G4int nHit; //number of times this pixel is hit
	G4int layerNo; //layer of photon origin

	//This flag is used to decide whether real hit is produced
	//It helps to discard background in simulation

public:
	void SetEvent(G4int n) 							{ event = n; }
	void SetPosition(const G4ThreeVector& pos) 		{ pixelPosition = pos; }
	void SetExactPosition(const G4ThreeVector& pos) { exactPosition = pos; }
	void SetPixelId(const G4int& id) 				{ pixelId  = id; }
	void SetVertexPos(const G4ThreeVector& pos) 	{ vertex = pos; }
	void SetDirection(const G4ThreeVector& dir) 	{ direction = dir; }
	void SetVertexDirection(const G4ThreeVector& dir) 	{ vertexDirection = dir; }
	void SetWL(G4double wl) 						{ wavelength = wl; }
	void Hit()			 							{ nHit++; }
	void SetLayerNo(const G4int i) 					{ layerNo = i; }

	G4int GetEvent() const						{ return event; }
	const G4ThreeVector& GetPosition() const	{ return pixelPosition; }
	const G4ThreeVector& GetExactPosition() const	{ return exactPosition; }
	const G4int& GetPixelId() const				{ return pixelId; }
	const G4ThreeVector& GetVertexPos() const	{ return vertex; }
	const G4ThreeVector& GetDirection() const 	{ return direction; }
	const G4ThreeVector& GetVertexDirection() const 	{ return vertexDirection; }
	G4double GetWL() const 						{ return wavelength; }
	G4bool IsHit() const 						{ return nHit?true:false; }
	G4int GetNhits() const 						{ return nHit; }
	G4int GetLayerNo() const					{ return layerNo; }
};

typedef G4THitsCollection<RichHit> RichHitsCollection;

extern G4Allocator<RichHit> RichHitAllocator;

inline void *RichHit::operator new(size_t)
{
	void *aHit;
	aHit = (void *) RichHitAllocator.MallocSingle();
	return aHit;
}

inline void RichHit::operator delete(void *aHit)
{
	RichHitAllocator.FreeSingle((RichHit *) aHit);
}

#endif
