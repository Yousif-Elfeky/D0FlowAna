/* **************************************************
 *
 *  Dzero pico ana.
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TVector3.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TNtuple.h"
#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "StMemStat.h"
#include "calmean.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#define piMass 0.13957
#define kMass  0.493677
#define massRes 0.013
#define nBuf 590
#define iBuf 2950 //nBuf * 5 events
#define nMaxBuffer 5
#define VzBinW 6.

//init cuts ------------------------------------------
#define yLocalCut    1.8   // |<1.8
//#define zLocalCut    3.2 // |<3.2

#define nsigtofshift 0.3
//#define kpiptcor     0.35  // kpt*pipt>0.35
//topo cuts
#define costCutOK    0     //false
#define cosstarCut   0.8
#define costPtCut    2.2
#ifndef DEBUG
#define DEBUG 1
#endif
ClassImp(StPicoDstarMixedMaker)
  StPicoDstarMixedMaker::StPicoDstarMixedMaker(char const * name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name), mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),mRunbyRunQVector(true)
{}

Int_t StPicoDstarMixedMaker::Init()
{
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  // -------------- USER VARIABLES -------------------------
    mFile = new TFile(mOutFileBaseName+".D0_v2.root", "RECREATE");
    fPion_up = new TF1("fPion_up", "(x<=1.6)*(5.1-2.25*x)+\(x>1.6)*1.5", 0.2, 20.);
    fPion_low = new TF1("fPion_low", "-2.", 0.2, 20.);
    fKaon_up = new TF1("fKaon_up", "(x<=2.5)*(6.129-1.9316*x)+\(x>2.5)*(1.3)", 0.2, 20.);
    fKaon_low = new TF1("fKaon_low", "(x<=1.7)*(-7.54+5.83*x-1.31*x*x)+\(x>1.7)*(-1.4149)", 0.2, 20.);

    for (int i = 0; i < nBuf; i++) {
        BufferEvent_NEvents[i] = 0;
        BufferEvent_Full[i] = 0;
    }
    
    for (int i = 0; i < iBuf; i++) {
        BufferEvent_nPions[i] = 0;
        BufferEvent_nKaons[i] = 0;
		BufferEvent_Id[i] = 0;
    }
    
    PI = 3.14159;
    twoPI = 6.28318;
    initHists();
    //openRecenterFile();
  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarMixedMaker::~StPicoDstarMixedMaker()
{}
//-----------------------------------------------------------------------------
void StPicoDstarMixedMaker::initHists(){
    int totalNum = 1590;
    if (mRunbyRunQVector){
        ifstream readnum;
        readnum.open(mRunNumList);
        int tmp=0;
        //if (DEBUG) cout<<"start initial run number..."<<endl;
        for (int i=0;i<totalNum;i++){
            readnum>>tmp;
            runnum.insert(pair<int,int>(tmp,i));
            
            //if (DEBUG) cout <<"run number : " << tmp <<" id :" <<runnum[tmp] <<endl;
        }
        readnum.close();
    }
  
    hnEvent = new TH1F("hnEvent", "# of events;Selection Flag", 10, 0., 10.);
    hVpdVz = new TH2F("hvpdVz", "vpdVz", 200, -100., 100., 200, -100., 100.);
    hVzDiff = new TH1F("hVzDiff", "", 200, -20., 20.);
    hRefMul = new TH1F("hRefMul", "hRefMul;RefMul", 700, 0., 700);
	hRefMul_weight = new TH1F("hRefMul_weight", "hRefMul_weight;hRefMul_weight", 700, 0., 700);
	hCentralityWeighted = new TH1F("hCentralityWeighted", "hCentralityWeighted", 9, 0, 9);
    hCentrality_noWgt = new TH1F("hCentrality_noWgt", "hCentrality_noWgt", 9, 0, 9);
    hnsigtofpi = new TH2F("hnsigtofpi", "nsigtofpi", 100, 0., 5., 2000, -10., 10.);
    hnsigtofk = new TH2F("hnsigtofk", "hnsigtofk", 100, 0., 5., 2000, -10., 10.);
    hEventPlane = new TH1F("hEventPlane", "hEventPlane", 170, -0.1, 3.3);
    
    /*hd0m = new TH1F("hd0m", "hd0m;d0m(GeV/c^{2})", 30000, 0., 3.);
    hd0m_like = new TH1F("hd0m_like", "hd0m_like;d0m(GeV/c^{2})", 30000, 0., 3.);
    hmix_d0m = new TH1F("hmix_d0m", "hmix_d0m;d0m(GeV/c^{2})", 30000, 0., 3.);
    hmix_d0m_like = new TH1F("hmix_d0m_like", "hmix_d0m_like;d0m(GeV/c^{2})", 30000, 0., 3.);
    hd0m_No = new TH1F("hd0m_No", "hd0m_No;d0m_No(GeV/c^{2})", 30000, 0., 3.);
    hd0m_like_No = new TH1F("hd0m_like_No", "hd0m_like_No;d0m_No(GeV/c^{2})", 30000, 0., 3.);
    hmix_d0m_No = new TH1F("hmix_d0m_No", "hmix_d0m_No;d0m_No(GeV/c^{2})", 30000, 0., 3.);
    hmix_d0m_like_No = new TH1F("hmix_d0m_like_No", "hmix_d0m_like_No;d0m_No(GeV/c^{2})", 30000, 0., 3.);
    */
    hmvsptcenPos = new TH3F("hmvsptcenPos", "m vs pt cen same event unlike", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hsmvsptcenPos = new TH3F("hsmvsptcenPos", "m vs pt cen same event like", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmmvsptcenPos = new TH3F("hmmvsptcenPos", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmsmvsptcenPos = new TH3F("hmsmvsptcenPos", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    /*hmvsptcenPos_No = new TH3F("hmvsptcenPos_No", "m vs pt cen same event unlike", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hsmvsptcenPos_No = new TH3F("hsmvsptcenPos_No", "m vs pt cen same event like", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmmvsptcenPos_No = new TH3F("hmmvsptcenPos_No", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmsmvsptcenPos_No = new TH3F("hmsmvsptcenPos_No", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);*/
    
    hmvsptcenNeg = new TH3F("hmvsptcenNeg", "m vs pt cen same event unlike", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hsmvsptcenNeg = new TH3F("hsmvsptcenNeg", "m vs pt cen same event like", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmmvsptcenNeg = new TH3F("hmmvsptcenNeg", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmsmvsptcenNeg = new TH3F("hmsmvsptcenNeg", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    /*hmvsptcenNeg_No = new TH3F("hmvsptcenNeg_No", "m vs pt cen same event unlike", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hsmvsptcenNeg_No = new TH3F("hsmvsptcenNeg_No", "m vs pt cen same event like", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmmvsptcenNeg_No = new TH3F("hmmvsptcenNeg_No", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);
    hmsmvsptcenNeg_No = new TH3F("hmsmvsptcenNeg_No", "m vs pt cen mix event", 9, 0, 9, 100, 0., 10., 300, 0., 3.);*/

}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Finish()
{
    mFile->cd();
    if (mRunbyRunQVector) {

	hnEvent->Write();
    hEventPlane->Write();
    hRefMul->Write();
	hRefMul_weight->Write();
    hCentralityWeighted->Write();
    hCentrality_noWgt->Write();
    hVpdVz->Write();
    hVzDiff->Write();
        
    hnsigtofpi->Write();
    hnsigtofk->Write();
        
    /*hd0m->Write();
    hmix_d0m->Write();
    hd0m_like->Write();
    hmix_d0m_like->Write();
    hd0m_No->Write();
    hmix_d0m_No->Write();
    hd0m_like_No->Write();
    hmix_d0m_like_No->Write();*/
        
    hmvsptcenPos->Write();
    hsmvsptcenPos->Write();
    hmmvsptcenPos->Write();
    hmsmvsptcenPos->Write();
    /*hmvsptcenPos_No->Write();
    hsmvsptcenPos_No->Write();
    hmmvsptcenPos_No->Write();
    hmsmvsptcenPos_No->Write();*/
    
    hmvsptcenNeg->Write();
    hsmvsptcenNeg->Write();
    hmmvsptcenNeg->Write();
    hmsmvsptcenNeg->Write();
    /*hmvsptcenNeg_No->Write();
    hsmvsptcenNeg_No->Write();
    hmmvsptcenNeg_No->Write();
    hmsmvsptcenNeg_No->Write();*/
        
  }
  mFile->Close();
    hnEvent->Delete();
    hRefMul->Delete();
	hRefMul_weight->Delete();
    hCentralityWeighted->Delete();
    hCentrality_noWgt->Delete();
    hVpdVz->Delete();
    hVzDiff->Delete();
    
    
    hnsigtofpi->Delete();
    hnsigtofk->Delete();
    
    
    /*hd0m->Delete();
    hmix_d0m->Delete();
    hd0m_like->Delete();
    hmix_d0m_like->Delete();
    hd0m_No->Delete();
    hmix_d0m_No->Delete();
    hd0m_like_No->Delete();
    hmix_d0m_like_No->Delete();*/
    
    hmvsptcenPos->Delete();
    hsmvsptcenPos->Delete();
    hmmvsptcenPos->Delete();
    hmsmvsptcenPos->Delete();
    /*hmvsptcenPos_No->Delete();
    hsmvsptcenPos_No->Delete();
    hmmvsptcenPos_No->Delete();
    hmsmvsptcenPos_No->Delete();*/
    
    hmvsptcenNeg->Delete();
    hsmvsptcenNeg->Delete();
    hmmvsptcenNeg->Delete();
    hmsmvsptcenNeg->Delete();
    /*hmvsptcenNeg_No->Delete();
    hsmvsptcenNeg_No->Delete();
    hmmvsptcenNeg_No->Delete();
    hmsmvsptcenNeg_No->Delete();*/
    hEventPlane->Delete();
    
    
  //if(fQProf) fQProf->Close();
  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Make()
{
  // StMemStat mem;
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  // -------------- USER ANALYSIS -------------------------
  StPicoEvent const * picoEvent = picoDst->event();
  hnEvent->Fill("all",1);
  mRunId = picoEvent->runId();
//   if (isBadrun(mRunId) || (!isGoodTrigger(picoEvent)))   return kStOK;
//   else 
  hnEvent->Fill("!Badrunlist & IsGoodTrigger",1);
  if (!isGoodEvent(picoEvent))  return kStOK;
  else hnEvent->Fill("isGoodEvent",1);
  float mVz = picoEvent->primaryVertex().z();
  float mVpdVz = picoEvent->vzVpd();
    
  //refmultCorrUtil = new StRefMultCorr("refmult","Isobar");
  StRefMultCorr *refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
  refmultCorrUtil -> init(mRunId);
  refmultCorrUtil->initEvent(picoEvent->refMult(),mVz,picoEvent->ZDCx());
  Bool_t isBadRun_Cen = refmultCorrUtil->isBadRun(mRunId);
  Bool_t isPileUpEvt_Cen = !refmultCorrUtil->passnTofMatchRefmultCut(picoEvent->refMult()*1.0,picoEvent->nBTOFMatch()*1.0);
  mCent  = refmultCorrUtil->getCentralityBin9();
  //if (isBadRun_Cen || isPileUpEvt_Cen || mCent<2 || mCent>8)  return kStOK;
  if (isBadRun_Cen || isPileUpEvt_Cen || mCent<0 || mCent>8)  return kStOK;
  else hnEvent->Fill("0-80%",1);
  CurrentEvent_Id = -999;
  CurrentEvent_Id = picoEvent -> eventId();
  
  weight = refmultCorrUtil->getWeight();
  mRefmult = picoEvent->refMult();
  hRefMul->Fill(mRefmult);
  hRefMul_weight->Fill(mRefmult, weight);
  hCentralityWeighted->Fill(mCent, weight);
  hCentrality_noWgt->Fill(mCent);
  hVpdVz->Fill(mVpdVz, mVz);
  hVzDiff->Fill(mVpdVz - mVz);
    
    //do bufferpointer
    if (mVz > -35. && mVz < -29.) CurrentEvent_Vz = 0;
    else if (mVz >= -29. && mVz < -23.) CurrentEvent_Vz = 1;
    else if (mVz >= -23. && mVz < -17.) CurrentEvent_Vz = 2;
    else if (mVz >= -17. && mVz < -11.) CurrentEvent_Vz = 3;
    else if (mVz >= -11. && mVz < -5.) CurrentEvent_Vz = 4;
    else if (mVz >= -5. && mVz < 1.) CurrentEvent_Vz = 5;
    else if (mVz >= 1. && mVz < 7.) CurrentEvent_Vz = 6;
    else if (mVz >= 7. && mVz < 13.) CurrentEvent_Vz = 7;
    else if (mVz >= 13. && mVz < 19.) CurrentEvent_Vz = 8;
    else                          CurrentEvent_Vz = 9;

	float eventPlane = calcEventPlane(picoDst, picoEvent, 2);
    hEventPlane->Fill(eventPlane);
	int mEventPlaneIndex = (int)(eventPlane / TMath::Pi() * 6);
    /*
    IMPORTANTE: The Buffer is a unique 3-digit number
    the fisrt digit from the left is the event plain index(0-5)
    the second digit is the centralty bin (0-8)
    the thrid digit is the vZ(0-9) each vZ is 6cm
    */
    BufferPointer = mEventPlaneIndex*100+mCent*10 + CurrentEvent_Vz; 
	CurrentEvent_nPions = 0;
    CurrentEvent_nKaons = 0;
    int nTracks = picoDst->numberOfTracks();
	for (int l = 0; l < nTracks; l++) {
        StPicoTrack* trk = picoDst->track(l);
        bool isprimary = trk -> isPrimary();
        if(!isprimary) continue;
        if (!isGoodTrack(trk, picoEvent,2)) continue;
        int index2tof = trk->bTofPidTraitsIndex();
        if(trk->isTofTrack()){
            StPicoBTofPidTraits const* const  tofPid = picoDst->btofPidTraits(index2tof);
            float tofylocal = tofPid->btofYLocal();
            if (fabs(tofylocal) > yLocalCut)  continue;
        }
        mcharge = trk->charge();
        nsigpi = trk->nSigmaPion();
        nsigk = trk->nSigmaKaon();
        
        TVector3 pmom = trk->pMom();
        TVector3 pV = picoEvent->primaryVertex();
        
        /*float pt = pmom.Perp();
        float eta = pmom.Eta();
        float phi = pmom.Phi();*/
        float pp = pmom.Mag();
        float beta = getTofBeta(trk);
        bool isTOFAvailable = (beta != std::numeric_limits<float>::quiet_NaN()) && beta > 0;
        if(pp<anaCuts::PIDPpCut && (!isTOFAvailable)) continue;
        nsigtofpi = (1. / beta - sqrt(piMass*piMass / pp / pp + 1)) / massRes;
        nsigtofk = (1. / beta - sqrt(kMass*kMass / pp / pp + 1)) / massRes;
        if (fabs(nsigpi) < 2.) {
            bool isPionCandidate = false;
            if (!isTOFAvailable)
                isPionCandidate = true;
            else {
                //hnsigtofpi->Fill(pp,nsigtofpi);
                if (nsigtofpi > fPion_low->Eval(pp) && nsigtofpi < fPion_up->Eval(pp)) isPionCandidate = true;
                else isPionCandidate = false;
            }
            
            if (isPionCandidate) {
                CurrentEvent_PionCharge[CurrentEvent_nPions] = mcharge;
                CurrentEvent_PionPx[CurrentEvent_nPions] = pmom.Px();
                CurrentEvent_PionPy[CurrentEvent_nPions] = pmom.Py();
                CurrentEvent_PionPz[CurrentEvent_nPions] = pmom.Pz();
                CurrentEvent_PionId[CurrentEvent_nPions] = l;
                CurrentEvent_nPions++;
            }
        }
        
        if (fabs(nsigk) < 2.) {
            bool isKaonCandidate = false;
            if (!isTOFAvailable) isKaonCandidate = true;
            else {
                //hnsigtofk->Fill(pp,nsigtofk);
                if (nsigtofk > fKaon_low->Eval(pp) && nsigtofk < fKaon_up->Eval(pp)) isKaonCandidate = true;
                else isKaonCandidate = false;
            }
            if (isKaonCandidate) {
                CurrentEvent_KaonCharge[CurrentEvent_nKaons] = mcharge;
                CurrentEvent_KaonPx[CurrentEvent_nKaons] = pmom.Px();
                CurrentEvent_KaonPy[CurrentEvent_nKaons] = pmom.Py();
                CurrentEvent_KaonPz[CurrentEvent_nKaons] = pmom.Pz();
                CurrentEvent_KaonId[CurrentEvent_nKaons] = l;
                CurrentEvent_nKaons++;
            }
        }
    }
    
	CurrentEvent_PairFound = false;
    makerealevent();
    if (CurrentEvent_PairFound) {
        makemixevent();
        copyToBuffer();
		hnEvent->Fill("PairFound",1);
    }
    
    //delete refmultCorrUtil;
    return kStOK;
  
}



bool StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
    for (auto trg : anaCuts::triggers)
    {
        if (picoEvent->isTrigger(trg)) return true;
    }
    
    return false;
}


bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, StPicoEvent const* picoEvent,int Index) const
{
    int index = Index;
    TVector3 pVtx = picoEvent->primaryVertex();
    TVector3 pmom = trk -> pMom();
    float pt = pmom.Perp();
    float eta = pmom.Eta();
    float nhits = trk->nHitsFit();
    int nhitmax = trk->nHitsMax();
    float ratio = fabs(nhits) / nhitmax;
    float dca = trk->gDCA(pVtx.x(), pVtx.y(), pVtx.z());
    if(index == 1){
        if (pt<=anaCuts::minPtCut1 || pt>=anaCuts::maxPtCut1) return false;
        else if (fabs(dca) > anaCuts::Dca1)   return false;
        else if (fabs(eta) >= anaCuts::EtaCut) return false;
        else if (fabs(nhits) < anaCuts::NHitsFit1) return false;
        else if (ratio<anaCuts::PtsRMin || ratio>anaCuts::PtsRMax) return false;
        else if (fabs(trk->charge())!=1) return false;
        else return true;
    }
    else {
        if (pt<=anaCuts::minPtCut2) return false;
        else if (fabs(dca) > anaCuts::Dca2)   return false;
        else if (fabs(eta) >= anaCuts::EtaCut) return false;
        else if (fabs(nhits) < anaCuts::NHitsFit2) return false;
        else if (ratio<anaCuts::PtsRMin || ratio>anaCuts::PtsRMax) return false;
        else if (fabs(trk->nHitsDedx()) <= anaCuts::mNHitsDedx) return false;
        else if (fabs(trk->charge())!=1) return false;
        else return true;
    }
}


bool StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const {
    
    float mRefmult = picoEvent->refMult();
    TVector3 pVtx = picoEvent->primaryVertex();
    if (pVtx.z() >= anaCuts::VzMax || pVtx.z() <= anaCuts::VzMin) return false;
    else if(fabs(pVtx.z() - picoEvent->vzVpd()) > anaCuts::VzDiff)  return false;
    else if (sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) > anaCuts::Vr) return false;
    else if (mRefmult < 0)  return false;
    else if ((fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror && fabs(pVtx.z()) < anaCuts::Verror)) return false;
    else return true;
}

float StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  float beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  } 
  return beta;
}

/*void StPicoDstarMixedMaker::openRecenterFile()
{
    fQProf = new TFile("StRoot/macros//recenter.histo.root");
    if(!fQProf->IsOpen())cout<<"###NOTE###: recenter file does not exist!!!"<<endl;
    prfQxCentEtaPluszplus = (TProfile2D*)fQProf->Get("prfQxCentEtaPluszplus");
    prfQyCentEtaPluszplus = (TProfile2D*)fQProf->Get("prfQyCentEtaPluszplus");
    prfQxCentEtaPluszminus = (TProfile2D*)fQProf->Get("prfQxCentEtaPluszminus");
    prfQyCentEtaPluszminus = (TProfile2D*)fQProf->Get("prfQyCentEtaPluszminus");
    
    prfQxCentEtaMinuszplus = (TProfile2D*)fQProf->Get("prfQxCentEtaMinuszplus");
    prfQyCentEtaMinuszplus = (TProfile2D*)fQProf->Get("prfQyCentEtaMinuszplus");
    prfQxCentEtaMinuszminus = (TProfile2D*)fQProf->Get("prfQxCentEtaMinuszminus");
    prfQyCentEtaMinuszminus = (TProfile2D*)fQProf->Get("prfQyCentEtaMinuszminus");

}
*/

void StPicoDstarMixedMaker::makerealevent()
{
    TLorentzVector kFourMom(0, 0, 0, 0);
    TLorentzVector piFourMom(0, 0, 0, 0);
    TLorentzVector d0FourMom(0, 0, 0, 0);
    for (int i = 0; i < CurrentEvent_nKaons; i++) {
        int   kaonid = CurrentEvent_KaonId[i];
        int   kaoncharge = CurrentEvent_KaonCharge[i];
        float kaonpx = CurrentEvent_KaonPx[i];
        float kaonpy = CurrentEvent_KaonPy[i];
        float kaonpz = CurrentEvent_KaonPz[i];
        kFourMom.SetXYZM(kaonpx, kaonpy, kaonpz, kMass);
        
        for (int j = 0; j < CurrentEvent_nPions; j++) {
            int pionid = CurrentEvent_PionId[j];
            if (kaonid == pionid) continue;
            int   pioncharge = CurrentEvent_PionCharge[j];
            float pionpx = CurrentEvent_PionPx[j];
            float pionpy = CurrentEvent_PionPy[j];
            float pionpz = CurrentEvent_PionPz[j];
            piFourMom.SetXYZM(pionpx, pionpy, pionpz, piMass);
            d0FourMom = kFourMom + piFourMom;
            float d0pt = d0FourMom.Pt();
            float d0m = d0FourMom.M();
            if(fabs(d0FourMom.Rapidity())>=anaCuts::D0yCut) continue;

            if ( kaoncharge * pioncharge < 0 ) {
                CurrentEvent_PairFound = true;
			//	hd0m->Fill(d0m,weight);
			//	hd0m_No->Fill(d0m);
            }

            //else {
              //  hd0m_like->Fill(d0m,weight);
               // hd0m_like_No->Fill(d0m);
           // }

            if (kaoncharge == -1 && pioncharge == 1) {
                    hmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                    //hmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                }

            else if (kaoncharge == 1 && pioncharge == -1) {
                    hmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                    //hmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                }

            else if (kaoncharge == 1 && pioncharge == 1) {
                    hsmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                    //hsmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                }

            else    {
                    hsmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                    //hsmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                }


        }//pi
    }//k
}

void StPicoDstarMixedMaker::makemixevent()
{
    TLorentzVector kFourMom(0, 0, 0, 0);
    TLorentzVector piFourMom(0, 0, 0, 0);
    TLorentzVector d0FourMom(0, 0, 0, 0);
    //int mCent = refmultCorrUtil->getCentralityBin9();
	for (int ic = 0; ic < CurrentEvent_nKaons; ic++) {
        int   kaoncharge = CurrentEvent_KaonCharge[ic];
        float kaonpx = CurrentEvent_KaonPx[ic];
        float kaonpy = CurrentEvent_KaonPy[ic];
        float kaonpz = CurrentEvent_KaonPz[ic];
        kFourMom.SetXYZM(kaonpx, kaonpy, kaonpz, kMass);
        
        for (int k = 0; k < BufferEvent_NEvents[BufferPointer]; k++) {
			if (BufferEvent_Id[k*nBuf + BufferPointer] == CurrentEvent_Id) continue;
            for (int ib = 0; ib < BufferEvent_nPions[k*nBuf + BufferPointer]; ib++) {
                int   pioncharge = BufferEvent_PionCharge[k*nBuf + BufferPointer][ib];
                float pionpx = BufferEvent_PionPx[k*nBuf + BufferPointer][ib];
                float pionpy = BufferEvent_PionPy[k*nBuf + BufferPointer][ib];
                float pionpz = BufferEvent_PionPz[k*nBuf + BufferPointer][ib];
                piFourMom.SetXYZM(pionpx, pionpy, pionpz, piMass);
                d0FourMom = kFourMom + piFourMom;
                float d0pt = d0FourMom.Pt();
                float d0m = d0FourMom.M();
                if(fabs(d0FourMom.Rapidity())>=anaCuts::D0yCut) continue;
                /*if( kaoncharge * pioncharge < 0) 
				{
					hmix_d0m->Fill(d0m,weight);
					hmix_d0m_No->Fill(d0m);
				}
                else  {
					hmix_d0m_like->Fill(d0m,weight);
					hmix_d0m_like_No->Fill(d0m);
				}*/

                if (kaoncharge == -1 && pioncharge == 1) {
                        hmmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                        //hmmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                    }
                else if (kaoncharge == 1 && pioncharge == -1) {
                        hmmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                        //hmmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                    }
                else if (kaoncharge == 1 && pioncharge == 1) {
                        hmsmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                        //hmsmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                    }
                else {
                        hmsmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                        //hmsmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                    }
            }  //Buffer Pion
        }//Buffer Event
    }//Current Kaon
    
    // current pi + buffer k
    for (int ic = 0; ic < CurrentEvent_nPions; ic++) {
        int   pioncharge = CurrentEvent_PionCharge[ic];
        float pionpx = CurrentEvent_PionPx[ic];
        float pionpy = CurrentEvent_PionPy[ic];
        float pionpz = CurrentEvent_PionPz[ic];
        piFourMom.SetXYZM(pionpx, pionpy, pionpz, piMass);
        for (int k = 0; k < BufferEvent_NEvents[BufferPointer]; k++) {
            if (BufferEvent_Id[k*nBuf + BufferPointer] == CurrentEvent_Id) continue;
            
            for (int ib = 0; ib < BufferEvent_nKaons[k*nBuf + BufferPointer]; ib++) {
                int kaoncharge = BufferEvent_KaonCharge[k*nBuf + BufferPointer][ib];
                float kaonpx = BufferEvent_KaonPx[k*nBuf + BufferPointer][ib];
                float kaonpy = BufferEvent_KaonPy[k*nBuf + BufferPointer][ib];
                float kaonpz = BufferEvent_KaonPz[k*nBuf + BufferPointer][ib];
                kFourMom.SetXYZM(kaonpx, kaonpy, kaonpz, kMass);
                d0FourMom = kFourMom + piFourMom;
                float d0pt = d0FourMom.Pt();
                float d0m = d0FourMom.M();
                if(fabs(d0FourMom.Rapidity())>=anaCuts::D0yCut) continue;
                
		/*		if ( kaoncharge * pioncharge < 0) {
					hmix_d0m->Fill(d0m,weight);
					hmix_d0m_No->Fill(d0m);
				}
                else  {
					hmix_d0m_like->Fill(d0m,weight);
					hmix_d0m_like_No->Fill(d0m);
				}*/

                if (kaoncharge == -1 && pioncharge == 1) {
                        hmmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                        //hmmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                    }
                else if (kaoncharge == 1 && pioncharge == -1) {
                        hmmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                        //hmmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                    }
                else if (kaoncharge == 1 && pioncharge == 1) {
                        hmsmvsptcenPos->Fill(mCent, d0pt, d0m,weight);
                        //hmsmvsptcenPos_No->Fill(mCent, d0pt, d0m);
                    }
                else {
                        hmsmvsptcenNeg->Fill(mCent, d0pt, d0m,weight);
                        //hmsmvsptcenNeg_No->Fill(mCent, d0pt, d0m);
                    }
            }//Buffer Kaon
        }//Buffer Event
    }//Current Pion
}

void StPicoDstarMixedMaker::copyToBuffer() {
	if (BufferEvent_NEvents[BufferPointer] >= nMaxBuffer) BufferEvent_Full[BufferPointer] = 1;
	TRandom3 *gRandom = new TRandom3(iran++);
    int eventPointer = -1;
	if (BufferEvent_Full[BufferPointer]) { // full - random rewrite one
        do {
            float rrr = gRandom->Rndm();
            eventPointer = (int)(nMaxBuffer*(1.0 - rrr));
        } while (eventPointer < 0 || eventPointer >= nMaxBuffer);
    }
    else { // not full
        eventPointer = BufferEvent_NEvents[BufferPointer];
    }
    delete gRandom;
    BufferEvent_Id[BufferPointer + nBuf * eventPointer] = CurrentEvent_Id;
    BufferEvent_nKaons[BufferPointer + nBuf * eventPointer] = CurrentEvent_nKaons;
    for (int ik = 0; ik < CurrentEvent_nKaons; ik++) {
        BufferEvent_KaonCharge[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_KaonCharge[ik];
        BufferEvent_KaonPx[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_KaonPx[ik];
        BufferEvent_KaonPy[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_KaonPy[ik];
        BufferEvent_KaonPz[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_KaonPz[ik];
    }
    BufferEvent_nPions[BufferPointer + nBuf * eventPointer] = CurrentEvent_nPions;
    for (int ik = 0; ik < CurrentEvent_nPions; ik++) {
        BufferEvent_PionCharge[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_PionCharge[ik];
        BufferEvent_PionPx[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_PionPx[ik];
        BufferEvent_PionPy[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_PionPy[ik];
        BufferEvent_PionPz[BufferPointer + nBuf * eventPointer][ik] = CurrentEvent_PionPz[ik];

    }
    
    if (BufferEvent_NEvents[BufferPointer] < nMaxBuffer) {
        BufferEvent_NEvents[BufferPointer]++;
    }

}


/*TVector2 StPicoDstarMixedMaker::QEtaGap(int iEta, int nEtaGaps) const  //calculate hadron v2, and D0 v2 need it
{
    TVector2 QEtaGap_(0, 0);
    int iEta_ = iEta;
    if (iEta_ < nEtaGaps) iEta_ = nEtaGaps - 1;
    if (iEta_ > 20 - nEtaGaps) iEta_ = 20 - nEtaGaps;
    for (int i = 0; i < 20; i++)
    {
        if (fabs(i - iEta_) >= nEtaGaps) QEtaGap_ += mQEta[i];
    }
    return QEtaGap_;
}*/

float StPicoDstarMixedMaker::calcEventPlane(StPicoDst const* const picoDst, StPicoEvent const* picoEvent, const int n) const
{
	float cossum_nocorrection = 0.;
	float sinsum_nocorrection = 0.;
	TVector3 pVtx = picoEvent->primaryVertex();
	const int nTrack = picoDst->numberOfTracks();
	for (int iTrack = 0; iTrack < nTrack; iTrack++) {
		StPicoTrack* mTrack = (StPicoTrack*)picoDst->track(iTrack);
		if (!mTrack) continue;
		const int nHitsFit = (int)fabs(mTrack->nHitsFit());
		if (nHitsFit < 20) continue;
		const float pt = mTrack->pMom().Perp();
		if (pt<0.2 || pt>2.0) continue; // pt<2 -- remove jet effect
		const float eta = mTrack->pMom().Eta();
		if (fabs(eta) > 1.0) continue;
		float dca = mTrack->gDCA(pVtx.x(), pVtx.y(), pVtx.z());
		if (fabs(dca) > 1.5) continue;
		const float phi = mTrack->pMom().Phi();
		const float cos_part_nocorrection = pt * cos(n*phi);
		const float sin_part_nocorrection = pt * sin(n*phi);
		cossum_nocorrection += cos_part_nocorrection;
		sinsum_nocorrection += sin_part_nocorrection;
	}
	TVector2 Q_nocorrection(cossum_nocorrection, sinsum_nocorrection);
	float eventPlane_nocorrection = Q_nocorrection.Phi();
	if (eventPlane_nocorrection < 0) eventPlane_nocorrection += 2.0*TMath::Pi();
	return eventPlane_nocorrection / n;

}
