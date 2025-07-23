#include <limits>
#include <cmath>

#include "StKaonPion.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

#include <iostream>

ClassImp(StKaonPion)


StKaonPion::StKaonPion(): mKaonIdx(std::numeric_limits<unsigned short>::min()),
						  mPionIdx(std::numeric_limits<unsigned short>::min()),
						  mPointingAngle(std::numeric_limits<float>::max()),
                          mDecayLength(std::numeric_limits<float>::min()),
                          mKaonDca(std::numeric_limits<float>::min()), 
                          mPionDca(std::numeric_limits<float>::min()),
                          mDcaDaughters(std::numeric_limits<float>::max()), 
                          mCosThetaStar(std::numeric_limits<float>::max()),
						  mMinv(std::numeric_limits<float>::min()),
						  mPt(std::numeric_limits<float>::min()),
						  mEta(std::numeric_limits<float>::min()),
						  mPhi(std::numeric_limits<float>::min())
{ }

StKaonPion::StKaonPion(StPicoTrack const& kaon, StPicoTrack const& pion, TVector3 const& vtx, float const bField) : StKaonPion()
{
   if (kaon.id() == pion.id())
   {
      return;
   }

   mKaonIdx = kaon.id(); // index of track in StPicoDstEvent
   mPionIdx = pion.id(); // index of track in StPicoDstEvent

   TLorentzVector mLorentzVector;

   /// prefixes code:
   ///   k means kaon
   ///   p means pion
   ///   kp means kaon-pion pair

   StPicoPhysicalHelix kHelix = kaon.helix(bField);
   StPicoPhysicalHelix pHelix = pion.helix(bField);

   // move origins of helices to the primary vertex origin
   kHelix.moveOrigin(kHelix.pathLength(vtx));
   pHelix.moveOrigin(pHelix.pathLength(vtx));

   // use straight lines approximation to get point of DCA of kaon-pion pair
   TVector3 const kMom = kHelix.momentum(bField * kilogauss);
   TVector3 const pMom = pHelix.momentum(bField * kilogauss);
   StPicoPhysicalHelix const kStraightLine(kMom, kHelix.origin(), 0, kaon.charge());
   StPicoPhysicalHelix const pStraightLine(pMom, pHelix.origin(), 0, pion.charge());
   
   pair<double, double> const ss = kStraightLine.pathLengths(pStraightLine);
   TVector3 const kAtDcaToPion = kStraightLine.at(ss.first);
   TVector3 const pAtDcaToKaon = pStraightLine.at(ss.second);

   // calculate DCA of pion to kaon at their DCA
   mDcaDaughters = (kAtDcaToPion - pAtDcaToKaon).Mag();

   // calculate Lorentz vector of kaon-pion pair
   TVector3 const kMomAtDca = kHelix.momentumAt(ss.first, bField * kilogauss);
   TVector3 const pMomAtDca = pHelix.momentumAt(ss.second, bField * kilogauss);

   TLorentzVector kFourMom;
   kFourMom.SetVectM(kMomAtDca,M_KAON_PLUS);
   TLorentzVector pFourMom;
   pFourMom.SetVectM(pMomAtDca,M_PION_PLUS);

   mLorentzVector = kFourMom + pFourMom;

   // calculate cosThetaStar
   TLorentzVector const kpFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
   TLorentzVector kFourMomStar(kFourMom);
   kFourMomStar.Boost(kpFourMomReverse.Vect());
   mCosThetaStar = std::cos(kFourMomStar.Vect().Angle(mLorentzVector.Vect()));

   // calculate pointing angle and decay length
   TVector3 const vtxToV0 = (kAtDcaToPion + pAtDcaToKaon) * 0.5 - vtx;
   mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
   mDecayLength = vtxToV0.Mag();

   // calculate DCA of tracks to primary vertex
   mKaonDca = (kHelix.origin() - vtx).Mag();
   mPionDca = (pHelix.origin() - vtx).Mag();

   TVector3 kMomPvtx  = kaon.gMom(vtx,bField);
   TVector3 piMomPvtx = pion.gMom(vtx,bField);

   mMinv = mLorentzVector.M();
   mPt 	 = mLorentzVector.Perp();
   mEta  = mLorentzVector.PseudoRapidity();
   mPhi  = mLorentzVector.Phi();

   mKaonPt  = kMomPvtx.Perp();
   mKaonEta = kMomPvtx.PseudoRapidity();
   mKaonPhi = kMomPvtx.Phi();
   mKaonDca = mKaonDca;
   mKaonCharge = kaon.charge();
   mKaonNHits = kaon.nHitsFit();
   mKaonNSigmaPion = kaon.nSigmaPion();
   mKaonNSigmaKaon = kaon.nSigmaKaon();
   mKaonDedx = kaon.dEdx();

   mcurvatureKa = kHelix.curvature();  // newly added
   mdipAngleKa = kHelix.dipAngle();
   mphaseKa = kHelix.phase();
   mhKa = kHelix.h();

   mPionPt  = piMomPvtx.Perp();
   mPionEta = piMomPvtx.PseudoRapidity();
   mPionPhi = piMomPvtx.Phi();
   mPionDca = mPionDca;
   mPionCharge = pion.charge();
   mPionNHits = pion.nHitsFit();
   mPionNSigmaPion = pion.nSigmaPion();
   mPionNSigmaKaon = pion.nSigmaKaon();
   mPionDedx = pion.dEdx();

   mcurvaturePi = pHelix.curvature();   // newly added
   mdipAnglePi = pHelix.dipAngle();
   mphasePi = pHelix.phase();
   mhPi = pHelix.h();
}
