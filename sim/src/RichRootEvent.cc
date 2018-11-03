#include <cstring>
#include "TClonesArray.h"
#include "RichRootEvent.hh"

ClassImp(RichRootParticle);
ClassImp(RichRootHit);
ClassImp(RichRootEvent);

TClonesArray* RichRootEvent::ghits = 0;

RichRootEvent::RichRootEvent() : TObject(),
	capacity(0),
	nhits(0),
	particle()
{
	if( !ghits ) ghits = new TClonesArray("RichRootHit",100);
    hits = ghits;
}

RichRootHit* RichRootEvent::AddHit()
{
	if( nhits>=(UInt_t)hits->GetSize() ) {
		hits->Expand(2*(nhits+1));
	}
	RichRootHit* p = new((*hits)[nhits++]) RichRootHit;
	return p;
}

