#ifndef StPicoD0Event__h
#define StPicoD0Event__h

/* **************************************************
 *  A specialized class for storing eventwise D0
 *  candidates. 
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  Modified by: Daniel Kikola
 *
 * **************************************************
 */

class StPicoEvent;

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"

class StKaonPion;
class StHaDron;

class StPicoD0Event : public TObject
{
public:
   StPicoD0Event();
   ~StPicoD0Event(){ clear("C");}
   void    clear(char const *option = "");
   void    addPicoEvent(StPicoEvent const& picoEvent);
   void    addKaonPion(StKaonPion const&);
    void    addHaDron(StHaDron const&); //Need to wrote StHaDron Class
   void    nKaons(int);
   void    nPions(int);
   void    refMultCorr(double);
   void    eventWeight(double);
   void    centrality(int);
    void    runId(int);
    void    eventId(int);
    void    bfield(float);
    void   vertex(TVector3);

   int   runId()   const;
   int   eventId() const;
    float  bfield() const;
   TClonesArray const* kaonPionArray()   const;
    TClonesArray const* HaDronArray()   const;
   int     nKaonPion()  const;
   int     nKaons() const;
   int     nPions() const;
    int     nHaDrons() const;
   double  refMultCorr();
   double  eventWeight();
   int	  centrality();
   TVector3 const& vertex() const;
    
private:
   // some variables below are kept in ROOT types to match the same ones in StPicoEvent
   int   mRunId;           // run number
   int   mEventId;         // event number
    float mBfield;


   double   mRefMultCorr;
   double   mEventWeight;
   int 	   mCentrality;
   TVector3 mVertex;
   int   mNKaonPion;       // number of stored pairs
    int   mNHaDron;       // number of stored HaDrons

   TClonesArray*        mKaonPionArray;
   static TClonesArray* fgKaonPionArray;
    
    TClonesArray*        mHaDronArray;
    static TClonesArray* fgHaDronArray;

   ClassDef(StPicoD0Event, 1)
};

inline void StPicoD0Event::refMultCorr(double n) { mRefMultCorr = n; }
inline void StPicoD0Event::eventWeight(double n) { mEventWeight = n; }
inline void StPicoD0Event::centrality(int n) { mCentrality = n; }
inline void StPicoD0Event::runId(int n) { mRunId = n; }
inline void StPicoD0Event::eventId(int n) { mEventId = n; }
inline void StPicoD0Event::bfield(float m) { mBfield = m; }
inline void StPicoD0Event::vertex(TVector3 t) { mVertex = t; }

inline TClonesArray const * StPicoD0Event::kaonPionArray()   const { return mKaonPionArray;}
inline TClonesArray const * StPicoD0Event::HaDronArray()   const { return mHaDronArray;}
inline int   StPicoD0Event::nKaonPion()  const { return mNKaonPion;}
inline int   StPicoD0Event::nHaDrons()  const { return mNHaDron;}
inline int StPicoD0Event::runId()   const { return mRunId; }
inline int StPicoD0Event::eventId() const { return mEventId; }
inline float StPicoD0Event::bfield() const { return mBfield; }
inline double StPicoD0Event::refMultCorr() { return mRefMultCorr; }
inline double StPicoD0Event::eventWeight() { return mEventWeight; }
inline int	StPicoD0Event::centrality() { return mCentrality; }

inline TVector3 const& StPicoD0Event::vertex() const { return mVertex; }
#endif
