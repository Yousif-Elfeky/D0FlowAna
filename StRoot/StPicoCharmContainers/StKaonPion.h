#ifndef StKaonPion_hh
#define StKaonPion_hh

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

class StKaonPion : public TObject
{
 public:
  StKaonPion();
  StKaonPion(StPicoTrack const& kaon, StPicoTrack const& pion, TVector3 const& vtx, float bField);
  ~StKaonPion() {}// please keep this non-virtual and NEVER inherit from this class 

//  StPicoPhysicalHelix kHelix;
//  StPicoPhysicalHelix pHelix;
  float m()    const;
  float pt()   const;
  float eta()  const;
  float phi()  const;
  float pointingAngle() const;
  float decayLength() const;
  float kaonDca() const;
  float pionDca() const;
  unsigned short   kaonIdx() const;
  unsigned short   pionIdx() const;
  float dcaDaughters() const;
  float cosThetaStar() const;
  float perpDcaToVtx() const;
  void  pionTofBeta(float);
  void  kaonTofBeta(float);
  unsigned short pionIdx();
  unsigned short kaonIdx();
  short pionCharge();
  short kaonCharge();
  short KaonNHits() const;
  short PionNHits() const;
    
    void  kaonOrigin(TVector3);
    void  pionOrigin(TVector3);
    TVector3 const& kaonOrigin() const; // used for access reduced data
    TVector3 const& pionOrigin() const; // used for access reduced data

   float KaonPt() const;
    float KaonEta() const;
    float KaonPhi() const;

    float PionPt() const;
    float PionEta() const;
    float PionPhi() const;


   void curvatureKaon(Double_t); // to save curvature values
   void dipAngleKaon(Double_t); // to save dip angle values
   void phaseKaon(Double_t); // to save phase values
   void hKaon(Int_t); // to save h parameter

   void curvaturePion(Double_t); // to save curvature values
   void dipAnglePion(Double_t); // to save dip angle values
   void phasePion(Double_t); // to save phase values
   void hPion(Int_t); // to save h parameter

  float curvatureKa();
  float dipAngleKa();
  float phaseKa();
  float hKa();

  float curvaturePi();
  float dipAnglePi();
  float phasePi();
  float hPi();



          
 private:
  unsigned short mKaonIdx; // index of track in StPicoDstEvent
  unsigned short mPionIdx; // index of track in StPicoDstEvent

  float mPointingAngle;
  float mDecayLength;
  float mKaonDca;
  float mPionDca;
  float mDcaDaughters;
  float mCosThetaStar;
  float mMinv;
  float mPt;
  float mEta;
  float mPhi;

  float mKaonPt;
  float mKaonEta;
  float mKaonPhi;
  short mKaonCharge;
  short mKaonNHits;
  float mKaonNSigmaPion;
  float mKaonNSigmaKaon;
  //short mKaonNSigmaProton;
  //short mKaonNSigmaElectron;
  float mKaonDedx;
  float mKaonTofBeta;

  float mPionPt;
  float mPionEta;
  float mPionPhi;
  short mPionCharge;
  short mPionNHits;
  float mPionNSigmaKaon;
  float mPionNSigmaPion;
  //short mPionNSigmaProton;
  //short mPionNSigmaElectron;
  float mPionDedx;
  float mPionTofBeta;
    
    TVector3 mKO;
    TVector3 mPO;

  Double_t mcurvatureKa;
  Double_t mdipAngleKa;
  Double_t mphaseKa;
  Int_t mhKa;

  Double_t mcurvaturePi;
  Double_t mdipAnglePi;
  Double_t mphasePi;
  Int_t mhPi;  


  ClassDef(StKaonPion,1)
};
inline float StKaonPion::m()    const { return mMinv; /*mLorentzVector.perp();*/}
inline float StKaonPion::pt()   const { return mPt;/*mLorentzVector.perp(); */ }
inline float StKaonPion::eta()  const { return mEta; /*mLorentzVector.pseudoRapidity(); */}
inline float StKaonPion::phi()  const { return mPhi; /* mLorentzVector.phi(); */}
inline float StKaonPion::pointingAngle() const { return mPointingAngle;}
inline float StKaonPion::decayLength() const { return mDecayLength;}
inline float StKaonPion::kaonDca() const { return mKaonDca;}
inline float StKaonPion::pionDca() const { return mPionDca;}
inline float StKaonPion::dcaDaughters() const { return mDcaDaughters;}
inline float StKaonPion::cosThetaStar() const { return mCosThetaStar;}
inline float StKaonPion::perpDcaToVtx() const { return mDecayLength*std::sin(mPointingAngle);}
inline void  StKaonPion::pionTofBeta(float n) { mPionTofBeta = n;}
inline void  StKaonPion::kaonTofBeta(float n) { mKaonTofBeta = n;}
inline unsigned short StKaonPion::pionIdx() { return mPionIdx; }
inline unsigned short StKaonPion::kaonIdx() { return mKaonIdx; }
inline short StKaonPion::pionCharge() { return mPionCharge; }
inline short StKaonPion::kaonCharge() { return mKaonCharge; }
inline void  StKaonPion::kaonOrigin(TVector3 t) { mKO = t;}
inline void  StKaonPion::pionOrigin(TVector3 t) { mPO = t;}

inline TVector3 const& StKaonPion::kaonOrigin() const { return mKO; }
inline TVector3 const& StKaonPion::pionOrigin() const { return mPO; }

inline float StKaonPion::KaonPt()  const { return mKaonPt; }
inline float StKaonPion::KaonEta() const { return mKaonEta; }
inline float StKaonPion::KaonPhi() const { return mKaonPhi; }
inline float StKaonPion::PionPt()  const { return mPionPt; }
inline float StKaonPion::PionEta() const { return mPionEta; }
inline float StKaonPion::PionPhi() const { return mPionPhi; }
inline short StKaonPion::KaonNHits() const {return mKaonNHits;}
inline short StKaonPion::PionNHits() const {return mPionNHits;}

inline void StKaonPion::curvatureKaon(Double_t a) {mcurvatureKa = a;}
inline void StKaonPion::dipAngleKaon(Double_t b) {mdipAngleKa = b;}
inline void StKaonPion::phaseKaon(Double_t c) {mphaseKa = c;}
inline void StKaonPion::hKaon(Int_t d)  {mhKa = d;}

inline void StKaonPion::curvaturePion(Double_t f) {mcurvaturePi = f;}
inline void StKaonPion::dipAnglePion(Double_t g) {mdipAnglePi = g;}
inline void StKaonPion::phasePion(Double_t k) {mphasePi = k;}
inline void StKaonPion::hPion(Int_t l)  {mhPi = l;}  
  
inline float StKaonPion::curvatureKa() {return mcurvatureKa;}
inline float StKaonPion::dipAngleKa() {return mdipAngleKa;}
inline float StKaonPion::phaseKa() {return mphaseKa;}
inline float StKaonPion::hKa()  {return mhKa;}

inline float StKaonPion::curvaturePi() {return mcurvaturePi;}
inline float StKaonPion::dipAnglePi() {return mdipAnglePi;}
inline float StKaonPion::phasePi() {return mphasePi;}
inline float StKaonPion::hPi()  {return mhPi;}

#endif

