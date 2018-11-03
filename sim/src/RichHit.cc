#include "RichHit.hh"

#include <G4ios.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4VisAttributes.hh>
#include <G4LogicalVolume.hh>
#include <G4Transform3D.hh>

#include "RichVisManager.hh"

G4Allocator <RichHit> RichHitAllocator;

RichHit::~RichHit() {}

void RichHit::Draw()
{
	G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
	if (pVVisManager && nHit) {
		G4Transform3D dummy;
		G4Circle circle(exactPosition);
		circle.SetScreenSize(3.0);
		circle.SetFillStyle(G4Circle::filled);
		G4VisAttributes attribs(G4Color(1.0,0.0,0.0));
		circle.SetVisAttributes(attribs);
		pVVisManager->Draw(circle,dummy);
	}
}
void RichHit::Print()
{
	G4cout<<"PMT hit: WL= "<<wavelength/nanometer<<"nm\n"
		  <<"         Position ("
		  <<pixelPosition.x()/cm<<"cm, "
		  <<pixelPosition.y()/cm<<"cm, "
		  <<pixelPosition.z()/cm<<"cm), nHit="<<nHit<<"\n";
}

// This is a forward declarations of an instantiated G4Allocator<Type> object.
// It has been added in order to make code portable for the GNU g++
// (release 2.7.2) compiler.
// Whenever a new Type is instantiated via G4Allocator, it has to be forward
// declared to make object code (compiled with GNU g++) link successfully.
//
#ifdef GNU_GCC
template class G4Allocator<RichHit>;
#endif
