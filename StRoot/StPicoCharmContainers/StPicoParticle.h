#ifndef StPicoParticle_hh
#define StPicoParticle_hh
#ifdef __ROOT__


#include <cmath>

#include "TObject.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StHaDron.h"
#include "StKaonPion.h"
using namespace std;
//class StHaDron
//class StKaonPion

class StPicoParticle : public StHaDron
{
public:
StPicoParticle(StPicoPhysicalHelix* helix, StHaDron* mHaParticle);
StPicoParticle(StPicoPhysicalHelix* helix, StKaonPion* mKPiPair);     
 
StPicoParticle(){}; // newly added
~StPicoParticle() {}

void CalculateTpcExitAndEntrancePoints(StPicoPhysicalHelix* helix);				  
int TpcLocalTransform(TVector3& aPoint, int& aSector, int& aRow, float& aU, double& aPhi);
double calcMergingPar(StPicoParticle* mTrcak1, StPicoParticle* mTrack2, double mMaxDuInner, double mMaxDuOuter, double mMaxDzInner, double mMaxDzOuter);
double NominalTpcAverageSeparation(StPicoParticle*, StPicoParticle* );

float ptKa() const;
float etaKa() const;
float phiKa() const;
unsigned short idxKa();
short chargeKa();

float ptPi() const;
float etaPi() const;
float phiPi() const;
unsigned short idxPi();
short chargePi();


void ptKa(float);
void etaKa(float);
void phiKa(float);
void idxKa(unsigned short);
void chargeKa(short);

void ptPi(float);
void etaPi(float);
void phiPi(float);
void idxPi(unsigned short);
void chargePi(short);

TVector3 PosSample(int i) const;        
//TVector3 PosSamplePlus(int j) const;
//TVector3 PosSampleMinus(int k) const;

private:
StPicoParticle* mTrack1;
StPicoParticle* mTrack2;
//StPicoParticle* MergPair;

float mPtKa;
float mEtaKa;
float mPhiKa;
unsigned short mIdxKa;
short mChargeKa;

float mPtPi;
float mEtaPi;
float mPhiPi;
unsigned short mIdxPi;
short mChargePi;



float mU[45];
float mZ[45];
int mSect[45];
TVector3 tmpTpcEntrancePoint; 
TVector3 tmpTpcExitPoint; 
TVector3 tmpPosSample[11];
//TVector3 tmpPosSamplePlus[11];
//TVector3 tmpPosSampleMinus[11];

/* inline float StPicoParticle::mU() {return mU[];}
inline float StPicoParticle::mZ() {return mZ[];}
inline float StPicoParticle::mSect() {return mSect[];} */
//inline float StPicoParticle::phi(float x) { return x;}

ClassDef(StPicoParticle,1)
};


inline float StPicoParticle::ptKa() const {return mPtKa;}
inline float StPicoParticle::etaKa() const {return mEtaKa;}
inline float StPicoParticle::phiKa() const {return mPhiKa;}
inline unsigned short StPicoParticle::idxKa() {return mIdxKa;}
inline short StPicoParticle::chargeKa() {return mChargeKa;}

inline float StPicoParticle::ptPi() const {return mPtPi;}
inline float StPicoParticle::etaPi() const {return mEtaPi;}
inline float StPicoParticle::phiPi() const {return mPhiPi;}
inline unsigned short StPicoParticle::idxPi() {return mIdxPi;}
inline short StPicoParticle::chargePi() {return mChargePi;}

inline void StPicoParticle::ptKa(float x) {mPtKa= x;}
inline void StPicoParticle::etaKa(float y) {mEtaKa = y;}
inline void StPicoParticle::phiKa(float z) {mPhiKa= z;}
inline void StPicoParticle::idxKa(unsigned short n) {mIdxKa= n;}
inline void StPicoParticle::chargeKa(short l) {mChargeKa= l;}

inline void StPicoParticle::ptPi(float a) {mPtPi= a;}
inline void StPicoParticle::etaPi(float b) {mEtaPi = b;}
inline void StPicoParticle::phiPi(float c) {mPhiPi= c;}
inline void StPicoParticle::idxPi(unsigned short d) {mIdxPi= d;}
inline void StPicoParticle::chargePi(short m) {mChargePi= m;}


inline TVector3 StPicoParticle::PosSample(int i) const {return tmpPosSample[i];}
//inline TVector3 StPicoParticle::PosSamplePlus(int j) const {return tmpPosSamplePlus[j];}
//inline TVector3 StPicoParticle::PosSampleMinus(int k) const {return tmpPosSampleMinus[k];}
//inline Double_t StPicoParticle::PosSample() const {return tmpPosSample[11];}
#endif
#endif
