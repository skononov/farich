/*
#include "RichUserTrackInformation.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichUserTrackInformation::RichUserTrackInformation() :
	status(active), nscatter(0), nrefl(0)
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
RichUserTrackInformation::~RichUserTrackInformation()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void RichUserTrackInformation::AddTrackStatusFlag(int s)
{
	if (s & active)				//track is now active
		status &= ~inactive;	//remove any flags indicating it is inactive
	else if (s & inactive)		//track is now inactive
		status &= ~active;		//remove any flags indicating it is active
	status |= s;				//add new flags
}
*/
