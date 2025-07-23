#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#include "StPicoD0Event.h"
#include "StKaonPion.h"
#include "StHaDron.h"

ClassImp(StPicoD0Event)

TClonesArray *StPicoD0Event::fgKaonPionArray = 0;
TClonesArray *StPicoD0Event::fgHaDronArray = 0;

//-----------------------------------------------------------------------
StPicoD0Event::StPicoD0Event() : mRunId(-1), mEventId(-1),
		mRefMultCorr (0), mEventWeight (-1), mCentrality (-1),
		mVertex(), mNKaonPion(0), mNHaDron(0), mKaonPionArray(NULL), mHaDronArray(NULL)
{
   if (!fgKaonPionArray) fgKaonPionArray = new TClonesArray("StKaonPion");
   mKaonPionArray = fgKaonPionArray;
    if (!fgHaDronArray) fgHaDronArray = new TClonesArray("StHaDron");
    mHaDronArray = fgHaDronArray;
}

//-----------------------------------------------------------------------
void StPicoD0Event::addPicoEvent(StPicoEvent const & picoEvent)
{
   // StPicoEvent variables
   mRunId = picoEvent.runId();
   mEventId = picoEvent.eventId();
   TVector3 primaryVtx = picoEvent.primaryVertex();

   mVertex.SetXYZ((Double_t)(primaryVtx.X()),(Double_t)(primaryVtx.Y()),(Double_t)(primaryVtx.Z()));
}

//-----------------------------------------------------------------------
void StPicoD0Event::clear(char const *option)
{
   mKaonPionArray->Clear(option);
    mHaDronArray->Clear(option);
   mRunId = -1;
   mEventId = -1;
    mBfield = -999.0;
   mVertex.SetXYZ(-999.,-999.,-999.);
   mNKaonPion = 0;
    mNHaDron = 0;
   mRefMultCorr = 0.0;
   mEventWeight = -1.0;
   mCentrality  = -1.0;


}
//---------------------------------------------------------------------
void StPicoD0Event::addKaonPion(StKaonPion const& t)
{
   new((*mKaonPionArray)[mNKaonPion++]) StKaonPion(t);
}

void StPicoD0Event::addHaDron(StHaDron const& v)
{
   new((*mHaDronArray)[mNHaDron++]) StHaDron(v);
}
