#include <limits>
#include <cmath>

#ifdef __ROOT__
#include "StHaDron.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

#include <iostream>

ClassImp(StHaDron)


StHaDron::StHaDron(): mPt(std::numeric_limits<float>::min()),
						  mEta(std::numeric_limits<float>::min()),
						  mPhi(std::numeric_limits<float>::min())
{ }

StHaDron::StHaDron(StPicoTrack const& trackin, TVector3 const& vrtx, float const bfield) : StHaDron()
{

   //TLorentzVector mLorentzVector;

   //mMinv = mLorentzVector.M();
   //mPt 	 = mLorentzVector.Perp();
   //mEta  = mLorentzVector.PseudoRapidity();
   //mPhi  = mLorentzVector.Phi();

/*    mPp = trackin.gMom(vrtx, bfield).Mag();
    mPt = trackin.gMom(vrtx, bfield).Perp();
    mEta = trackin.gMom(vrtx, bfield).PseudoRapidity();
    mPhi = trackin.gMom(vrtx, bfield).Phi();      */

     mPp = trackin.pMom().Mag();
    mPt = trackin.pMom().Perp();
    mEta = trackin.pMom().PseudoRapidity();
    mPhi = trackin.pMom().Phi();

    mCharge = trackin.charge();
    mNSigmaPion = trackin.nSigmaPion();
    mNSigmaKaon = trackin.nSigmaKaon();
    mNSigmaProton = trackin.nSigmaProton();
    mNSigmaElectron = trackin.nSigmaElectron();
    mDedx = trackin.dEdx();
    mNHits = trackin.nHitsFit();
    mNFitPoss = trackin.nHitsPoss();
    mIdx = trackin.id();
}
#endif // __ROOT__
