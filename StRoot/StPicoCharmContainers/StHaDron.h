#ifndef StHaDron_hh
#define StHaDron_hh
#ifdef __ROOT__

/* **************************************************
 *  A specialized pair class for calculating K-π pair 
 *  lorentz vector and topological decay parameters 
 *  and storing them.
 *
 *  Authors:  Xin Dong (xdong@lbl.gov),
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  Modified by D. Kikola;
 *
 * **************************************************
 */
#include <cmath>

#include "TObject.h"
#include "TVector3.h"
#include "TClonesArray.h"

class StPicoTrack;
class StPicoEvent;

class StHaDron : public TObject
{
 public:
  StHaDron();
  //StHaDron(StPicoTrack const& kaon, StPicoTrack const& pion, TVector3 const& vtx, float bField);
    StHaDron(StPicoTrack const& trackin, TVector3 const& vrtx, float const bfield);
  ~StHaDron() {}// please keep this non-virtual and NEVER inherit from this class

  //float px()    const;
   // float py()    const;
    float pp()    const;
    float pt()   const;
    float eta()  const;
    float phi()  const;
    short charge();
    short nhits() const;
    short nfitposs() const;
    float sigmaPion() const;
    float sigmaKaon() const;
    float sigmaProton() const;
    float sigmaElectron() const;
    float dedx() const;
    float dca() const;      // used for D0 reduced tree production
    float tofBeta() const;
    void pt(float); // used to set pt values
    void eta(float); // used to set eta values
    void phi(float); // used to set phi values
    void idx(unsigned short); // used to set idx values
    void chargeHa(short); // used to set charge values
    void pp(float); // used to store momentum values
    
    unsigned short Idx();
    
    void  tofBeta(float);
    void  dca(double);     // used for access reduced data
    
    void  hadronOrigin(TVector3);   // used for D0 reduced tree production
    TVector3 const& hadronOrigin() const; // used for access reduced data

   void curvature(Double_t); // to save curvature values
   void dipAngle(Double_t); // to save dip angle values
   void phase(Double_t); // to save phase values
   void h(Int_t); // to save h parameter

   float curvature();
   float dipAngle();
   float phase();
   int h();

    
 private:
    
    unsigned short mIdx; // index of track in StPicoDstEvent

    float mPp;
    float mPt;
    float mEta;
    float mPhi;
    short mCharge;
    short mNHits;
    short mNFitPoss;
    float mNSigmaPion;
    float mNSigmaKaon;
    float mNSigmaProton;
    float mNSigmaElectron;
    float mDedx;
    float mTofBeta;
    double mDca;
    TVector3 mHO;
   
    Double_t mcurvature;
    Double_t mdipAngle;
    Double_t mphase;
    Int_t mh;  

    
  ClassDef(StHaDron,1)
};
//inline float StHaDron::px()    const { return mPx; /*mLorentzVector.perp();*/}
//inline float StHaDron::py()    const { return mPy; /*mLorentzVector.perp();*/}
inline float StHaDron::pp()    const { return mPp; /*mLorentzVector.perp();*/}
inline float StHaDron::pt()   const { return mPt;/*mLorentzVector.perp(); */ }
inline float StHaDron::eta()  const { return mEta; /*mLorentzVector.pseudoRapidity(); */}
inline float StHaDron::phi()  const { return mPhi; /* mLorentzVector.phi(); */}
inline void  StHaDron::tofBeta(float n) { mTofBeta = n;}
//inline unsigned short StHaDron::Idx() { return mIdx; }
inline short StHaDron::charge() { return mCharge; }
inline float StHaDron::sigmaPion()   const { return mNSigmaPion;/*mLorentzVector.perp(); */ }
inline float StHaDron::sigmaKaon()   const { return mNSigmaKaon;/*mLorentzVector.perp(); */ }
inline float StHaDron::sigmaProton()   const { return mNSigmaProton;/*mLorentzVector.perp(); */ }
inline float StHaDron::sigmaElectron()   const { return mNSigmaElectron;/*mLorentzVector.perp(); */ }
inline float StHaDron::dedx()   const { return mDedx;/*mLorentzVector.perp(); */ }
inline short StHaDron::nhits()   const { return mNHits;/*mLorentzVector.perp(); */ }
inline short StHaDron::nfitposs()   const { return mNFitPoss;/*mLorentzVector.perp(); */ }
inline void  StHaDron::dca(double m) { mDca = m;}
inline unsigned short StHaDron::Idx() { return mIdx; }
inline void  StHaDron::hadronOrigin(TVector3 t) { mHO = t;}
inline void StHaDron::pt(float x) {mPt =x;}
inline void StHaDron::eta(float y) {mEta = y;}
inline void StHaDron::phi(float z) {mPhi= z;}
inline void StHaDron::idx(unsigned short a) {mIdx= a;}
inline void StHaDron::chargeHa(short k) {mCharge= k;}
inline void StHaDron::pp(float l) {mPp= l;}

inline float StHaDron::tofBeta() const {return mTofBeta;}
inline float StHaDron::dca() const {return mDca; }
inline TVector3 const& StHaDron::hadronOrigin() const { return mHO; }

inline void StHaDron::curvature(Double_t a) {mcurvature = a;}
inline void StHaDron::dipAngle(Double_t b) {mdipAngle = b;}
inline void StHaDron::phase(Double_t c) {mphase = c;}
inline void StHaDron::h(Int_t d)  {mh = d;}  

inline float StHaDron::curvature() {return mcurvature;}
inline float StHaDron::dipAngle() {return mdipAngle;}
inline float StHaDron::phase() {return mphase;}
inline int StHaDron::h() {return mh;}

#endif
#endif

