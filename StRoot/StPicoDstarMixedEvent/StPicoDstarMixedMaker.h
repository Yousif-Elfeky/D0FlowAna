#ifndef StPicoDstarMixedMaker_h
#define StPicoDstarMixedMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *  This can also be extended to analyze Dzero v2.
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h" // --- ADDED FOR V_N ANALYSIS ---
#include "TVector2.h" // --- ADDED FOR V_N ANALYSIS ---
#include "TF1.h"
#include "TVector3.h"
class TString;
class TFile;
class TNtuple;
class TProfile;
class TProfile2D;
class TProfile3D;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;
class StPicoDst;
class StRefMultCorr;
class CentralityMaker;

class StPicoDstarMixedMaker : public StMaker
{
  public:
    StPicoDstarMixedMaker(char const * name, TString const inputFilesList,
                          TString const outBaseName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoDstarMixedMaker();
    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
  
  private:
    StPicoDstarMixedMaker() {}
    void initHists();
    void  openRecenterFile ( );
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodEvent(StPicoEvent const* const picoEvent) const;
    bool isGoodTrack(StPicoTrack const* trk,StPicoEvent const* picoEvent,int Index) const;
    float getTofBeta(StPicoTrack const* const trk) const;
    void makerealevent();
    void makemixevent();
    void copyToBuffer();
    float calcEventPlane(StPicoDst const* const picoDst, StPicoEvent const* picoEvent, const int n) const;
	  void calculateHadronV2(StPicoEvent const* const picoEvent,std::map<int,int> runnum) const;
    TVector2 QEtaGap(int iEta, int nEtaGaps) const;
    StPicoDstMaker* mPicoDstMaker;
    StRefMultCorr* refmultCorrUtil;
    TString mInputFilesList;
    TString mOutFileBaseName;
    bool isBadrun(Int_t runId);
    // --- ADDED FOR V_N ANALYSIS ---
    void getQVectors(StPicoDst const* picoDst, TVector2 Q[3], int n) const;
    void calculateEventPlaneResolution();
    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate 

  public:
    void setRunNumList(string list){
      mRunNumList = list;
    }
    void setRunbyRunQVector(bool b){
      mRunbyRunQVector = b;
    }

    void getBadruns(string inputFileName);
  private: 
  
    TFile* mFile;

    std::map<int,int> runnum;
    string mRunNumList;
    vector<int> mBadRun;
    
    TH1F*      hCentralityWeighted;
    TH1F*      hCentrality_noWgt;
    //TProfile*  prfCosResolutionEtaCent;
    TH2F*         hVpdVz;
    TH2F*         hnsigtofpi;
    TH2F*         hnsigtofk;
    
    TH1F* hnEvent;
    TH1F* hEventPlane;
    TH1F* hRefMul;
    TH1F* hRefMul_weight;
    TH1F* hVzDiff;
    TH1F* hd0m;
    TH1F* hd0m_like;
    TH1F* hmix_d0m;
    TH1F* hmix_d0m_like;
    TH1F* hd0m_No;
    TH1F* hd0m_like_No;
    TH1F* hmix_d0m_No;
    TH1F* hmix_d0m_like_No;
    
    
    TH3F* hmvsptcenPos;
    TH3F* hsmvsptcenPos;
    TH3F* hmmvsptcenPos;
    TH3F* hmsmvsptcenPos;
    TH3F* hmvsptcenPos_No;
    TH3F* hsmvsptcenPos_No;
    TH3F* hmmvsptcenPos_No;
    TH3F* hmsmvsptcenPos_No;
    
    TH3F* hmvsptcenNeg;
    TH3F* hsmvsptcenNeg;
    TH3F* hmmvsptcenNeg;
    TH3F* hmsmvsptcenNeg;
    TH3F* hmvsptcenNeg_No;
    TH3F* hsmvsptcenNeg_No;
    TH3F* hmmvsptcenNeg_No;
    TH3F* hmsmvsptcenNeg_No;

    //int byDayOrRun;    //0: by run;  1: by day
    //int timePointer;   //
    int  mRunId;
    int  mCent;
    float weight;
    float mRefmult;

    // --- ADDED FOR V_N ANALYSIS ---
    TVector2 mEventPlaneV1[3]; // 0: eta<neg_gap, 1: eta>pos_gap, 2: full
    TVector2 mEventPlaneV2[3];
    TProfile *hCos_v1_ab, *hCos_v1_ac, *hCos_v1_bc;
    TProfile *hCos_v2_ab, *hCos_v2_ac, *hCos_v2_bc;
    TH1F *hEventPlaneRes_v1, *hEventPlaneRes_v2;
    THnSparseF *hD0_v1_UL, *hD0_v1_LS, *hD0_v1_MixUL, *hD0_v1_MixLS;
    THnSparseF *hD0_v2_UL, *hD0_v2_LS, *hD0_v2_MixUL, *hD0_v2_MixLS;

    //float M_pion=0.139570;//GeV
    float PI;
    float twoPI;
    bool mRunbyRunQVector;
    bool mDaybyDayQVector;
     
    TFile*        fQProf;
    
    //D0_v2
    TF1* fPion_low;
    TF1* fPion_up;
    TF1* fKaon_low;
    TF1* fKaon_up;
    
    int mcharge;
    float nsigpi;
    float nsigk;
    float nsigtofpi;
    float nsigtofk;
    int   mEventPlaneIndex;
    
    int   CurrentEvent_nPions;               //Pion signal in CurrentEvent
    int   CurrentEvent_PionCharge[1000];
    float CurrentEvent_PionPx[1000];
    float CurrentEvent_PionPy[1000];
    float CurrentEvent_PionPz[1000];
    int   CurrentEvent_PionId[1000];
    int   CurrentEvent_nKaons;               //Kaon signal
    int   CurrentEvent_KaonCharge[1000];
    float CurrentEvent_KaonPx[1000];
    float CurrentEvent_KaonPy[1000];
    float CurrentEvent_KaonPz[1000];
    int   CurrentEvent_KaonId[1000];
    
    bool  CurrentEvent_PairFound;
    
    int   CurrentEvent_Vz;
    int   BufferPointer;
    int   CurrentEvent_Id;
    int   BufferEvent_NEvents[590];
    int   iran;
    int   BufferEvent_Full[590];
    int   BufferEvent_Id[2950];
    int   BufferEvent_nPions[2950];
    int   BufferEvent_nKaons[2950];
    //float BufferEvent_Qx[2950];
    //float BufferEvent_Qy[2950];

    int   BufferEvent_PionCharge[2950][1000];
    float BufferEvent_PionPx[2950][1000];
    float BufferEvent_PionPy[2950][1000];
    float BufferEvent_PionPz[2950][1000];
    
    int   BufferEvent_KaonCharge[2950][1000];
    float BufferEvent_KaonPx[2950][1000];
    float BufferEvent_KaonPy[2950][1000];
    float BufferEvent_KaonPz[2950][1000];

    bool mRunbyRunQA;
        ClassDef(StPicoDstarMixedMaker, 1)
};

inline void StPicoDstarMixedMaker::getBadruns(string inputFileName){
    ifstream fin(inputFileName.c_str());
    if(!fin){
      cout <<"no Bad runs list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId ) {
      mBadRun.push_back(runId);
    }
    cout << "get Bad runs list [OK]" << endl;
}
inline  bool StPicoDstarMixedMaker::isBadrun(Int_t runId){
    vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), runId);
    return ( iter != mBadRun.end() ) ;
}

#endif
