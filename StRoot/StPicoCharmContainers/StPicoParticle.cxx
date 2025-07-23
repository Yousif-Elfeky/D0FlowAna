#include <limits>
#include <cmath>

#ifdef __ROOT__

#include "StPicoParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"


#include <iostream>

using namespace std;

ClassImp(StPicoParticle) 

StPicoParticle::StPicoParticle(StPicoPhysicalHelix* helix, StHaDron* mHaParticle){

	CalculateTpcExitAndEntrancePoints(helix);
	pt(mHaParticle->pt());
	eta(mHaParticle->eta());
	phi(mHaParticle->phi());
	idx(mHaParticle->Idx());
	chargeHa(mHaParticle->charge());
	pp(mHaParticle->pp());
			    
}

StPicoParticle::StPicoParticle(StPicoPhysicalHelix* helix, StKaonPion* mKPiPair){

	CalculateTpcExitAndEntrancePoints(helix);
	ptPi(mKPiPair->PionPt());
	ptKa(mKPiPair->KaonPt());
	etaPi(mKPiPair->PionEta());
        etaKa(mKPiPair->KaonEta());
	phiPi(mKPiPair->PionPhi());
        phiKa(mKPiPair->KaonPhi());
	idxPi(mKPiPair->pionIdx());
	idxKa(mKPiPair->kaonIdx());
	chargePi(mKPiPair->pionCharge());
	chargeKa(mKPiPair->kaonCharge());
	
}

/* StPicoParticle::StPicoParticle(StPicoPhysicalHelix* helix, StKaonPion* mKa){

        CalculateTpcExitAndEntrancePoints(helix);
        pt(mKPiPair->KaonPt());
//      pt(mKPiPair->pt());
      //        eta(mKPiPair->eta());
              eta(mKPiPair->KaonEta());
                      phi(mKPiPair->KaonPhi());
                      //        phi(mKPiPair->phi());
                      //      idx(mKPiPair->kaonIdx());
                              idx(mKPiPair->kaonIdx());
                              //      pt(mKPiPair->KaonPt());
                              //      pt(mKPiPair->PionPt());

                              }*/


// calculate TPC exit and entrance points
void StPicoParticle::CalculateTpcExitAndEntrancePoints(StPicoPhysicalHelix* helix){
 //template< class T1, class T2 > struct pair;
 
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */
//StHaDron* mHa;
//cout << " phi angle of hadron " << mHa->phi() << endl;
 std::pair<double, double> candidates;
 
  double sideLength;  // this is how much length to go to leave through sides of TPC
  double endLength;  // this is how much length to go to leave through endcap of TPC
  // figure out how far to go to leave through side...
  candidates = helix->pathLength(200.0);  // bugfix MAL jul00 - 200cm NOT 2cm // pathlength at radial distance r
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second; 

  static TVector3 WestEnd(0.,0.,200.);  // bugfix MAL jul00 - 200cm NOT 2cm
  static TVector3 EastEnd(0.,0.,-200.); // bugfix MAL jul00 - 200cm NOT 2cm
  static TVector3 EndCapNormal(0.,0.,1.0);

  endLength = helix->pathLength(WestEnd,EndCapNormal);  // pathlength at intersection with plane
  if (endLength < 0.0) endLength = helix->pathLength(EastEnd,EndCapNormal);   // pathlength at intersection with plane

  if (endLength < 0.0) cout << 
			 "D0AzimuthalCorAnalyzer::CalculateTpcExitAndEntrancePoints(): "
                            << "Hey -- I cannot find an exit point out endcaps" << endl;
  // OK, firstExitLength will be the shortest way out of the detector...
  double firstExitLength = (endLength < sideLength) ? endLength : sideLength;
  // now then, let's return the POSITION at which particle leaves TPC...
   tmpTpcExitPoint = helix->at(firstExitLength);
  // Finally, calculate the position at which the track crosses the inner field cage
  
   // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */

  candidates = helix->pathLength(50.0);  // bugfix MAL jul00 - 200cm NOT 2cm
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
  
  //  cout << "sideLength 2 ="<<sideLength << endl;
   tmpTpcEntrancePoint = helix->at(sideLength);
  // This is the secure way !  
  if (std::isnan(tmpTpcEntrancePoint.x()) || 
      std::isnan(tmpTpcEntrancePoint.y()) || 
      std::isnan(tmpTpcEntrancePoint.z()) ){ 
    cout << "tmpTpcEntrancePoint NAN"<< endl; 
   // cout << "tmpNominalTpcEntrancePoint = " <<tmpTpcEntrancePoint<< endl;
    tmpTpcEntrancePoint.SetX(-9999.);
    tmpTpcEntrancePoint.SetY(-9999.);
    tmpTpcEntrancePoint.SetZ(-9999.);
  } 
  	
  if (std::isnan(tmpTpcExitPoint.x()) || 
      std::isnan(tmpTpcExitPoint.y()) || 
      std::isnan(tmpTpcExitPoint.z()) ) {
//     cout << "tmpTpcExitPoint NAN set at (-9999,-9999,-9999)"<< endl; 
//     cout << "tmpTpcExitPoint X= " <<tmpTpcExitPoint->x()<< endl;
//     cout << "tmpTpcExitPoint Y= " <<tmpTpcExitPoint->y()<< endl;
//     cout << "tmpTpcExitPoint Z= " <<tmpTpcExitPoint->z()<< endl;
    tmpTpcExitPoint.SetX(-9999.);
    tmpTpcExitPoint.SetY(-9999.);
    tmpTpcExitPoint.SetZ(-9999.);
  }
  
  // calculate the "nominal" position at N radii (say N=11) 
  // within the TPC, and for a pair cut
  // use the average separation of these N
  int irad = 0;
  
   // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */

  candidates = helix->pathLength(50.0);
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
  while (irad<11 && !std::isnan(sideLength)){
    float radius = 50.0 + irad*15.0;
    
     // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */

    candidates = helix->pathLength(radius);
    sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    tmpPosSample[irad] = helix->at(sideLength);
   // Double_t tmpPosSample [irad] = (Double_t)tmpPosSample;
    if(std::isnan(tmpPosSample[irad].x()) ||
       std::isnan(tmpPosSample[irad].y()) ||
       std::isnan(tmpPosSample[irad].z()) 
       ){
      cout << "tmpPosSample for radius=" << radius << " NAN"<< endl; 
     // cout << "tmpPosSample=(" << tmpPosSample[irad] << ")"<< endl;
      tmpPosSample[irad] =  TVector3(-9999.,-9999.,-9999);
    }
    irad++;
    if (irad<11){
      float radius = 50.0 + irad*15.0;
      
       // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */

      candidates = helix->pathLength(radius);  // pathlength at cylindrical r
      sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    }
   }
   for (int i = irad; i<11; i++)
     {
       tmpPosSample[i] =  TVector3(-9999.,-9999.,-9999);   
     }
     
     static float tRowRadius[45] = {60,64.8,69.6,74.4,79.2,84,88.8,93.6,98.8, 
				 104,109.2,114.4,119.6,127.195,129.195,131.195,
				 133.195,135.195,137.195,139.195,141.195,
				 143.195,145.195,147.195,149.195,151.195,
				 153.195,155.195,157.195,159.195,161.195,
				 163.195,165.195,167.195,169.195,171.195,
				 173.195,175.195,177.195,179.195,181.195,
				 183.195,185.195,187.195,189.195};
  int tRow,tSect,tOutOfBound;
  double tLength,tPhi;
  float tU;
  TVector3 tPoint;
  TVector3 tn(0,0,0);
  TVector3 tr(0,0,0);
  int ti =0;
  /*float mU[45];
  float mZ[45];
  int mSect[45];*/

  // test to enter the loop
   // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl;*/

  candidates =  helix->pathLength(tRowRadius[ti]);
  tLength = (candidates.first > 0) ? candidates.first : candidates.second;
  if (std::isnan(tLength)){
    cout <<"tLength Init tmp NAN" << endl;
   // cout <<"padrow number= "<<ti << "not reached" << endl;
    cout << "*** DO NOT ENTER THE LOOP***" << endl;
    mSect[ti]=-1;//sector
  }
  // end test
  while(ti<45 && !std::isnan(tLength)){
  // checking helix parameters
/* cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl; */

    candidates =  helix->pathLength(tRowRadius[ti]);
    //cout << "padrow raduius at " << ti << " is " << " " << tRowRadius[ti] << endl; // NEWLY ADDED
   // cout << " candidates at row " << ti << " are " << candidates.first << " and " << candidates.second << endl;
    tLength = (candidates.first > 0) ? candidates.first : candidates.second;
  //  cout << "tLength at " << ti << " is " << tLength << endl;
    
    if (std::isnan(tLength)){
      cout <<"tLength loop 1st NAN" << endl;
    //  cout <<"padrow number=  " << ti << " not reached" << endl;
      cout << "*** THIS IS AN ERROR SHOULDN'T  LOOP ***" << endl;
      mSect[ti]=-1;//sector
    }
    tPoint = helix->at(tLength);
    //cout <<"hit points are= "<< tPoint << endl; // NEW LINE
    // Find which sector it is on
      
  // cout << " global x, y, z coordinates are " << tPoint.x() << " " << tPoint.y() << " " << tPoint.z() << " at sector" << ti << endl;
    TpcLocalTransform(tPoint, mSect[ti], tRow, tU, tPhi);
    if (std::isnan(mSect[ti])){
      cout <<"***ERROR mSect"<< endl; 
    }
    if (std::isnan(tRow)){
      cout <<"***ERROR tRow"<< endl;
    }
    if (std::isnan(tU)){
      cout <<"***ERROR tU"<< endl;
    }
    if (std::isnan(tPhi)){
      cout <<"***ERROR tPhi"<< endl;
    }  
    // calculate crossing plane
    tn.SetX(cos(tPhi));
    tn.SetY(sin(tPhi));       
    tr.SetX(tRowRadius[ti]*cos(tPhi));
    tr.SetY(tRowRadius[ti]*sin(tPhi));
    // find crossing point
    tLength = helix->pathLength(tr,tn); // This pathlength defined the pathlength at intersection with plane
    if (std::isnan(tLength)){
      cout <<"tLength loop 2nd  NAN" << endl;
      cout <<"padrow number=  " << ti << " not reached" << endl;
      mSect[ti]=-2;//sector
    }
    tPoint = helix->at(tLength);
    mZ[ti] = tPoint.z();
  //  cout <<"coordinates along Z axis are= " << mZ[ti] << endl; // NEW LINE
    tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow, mU[ti],tPhi);
   // vector <vector <vector <TVector3 > > >tPointVec; // NEWLY ADDED >>> Vector defination to store tPoint
   // tPointVec.push_back(tPoint);  // NEWLY ADDED >>> info storing 
    if (std::isnan(tSect)){
      cout <<"***ERROR tSect 2"<< endl; 
    }
    if (std::isnan(tRow)){
      cout <<"***ERROR tRow 2"<< endl;
    }
    if (std::isnan(mU[ti])){
      cout <<"***ERROR mU[ti] 2"<< endl;
    }
    if (std::isnan(tPhi)){
      cout <<"***ERROR tPhi 2 "<< endl;
    }  
    if(tOutOfBound || (mSect[ti] == tSect && tRow!=(ti+1))){
      mSect[ti]=-2;
      //	  cout << "missed once"<< endl;
    }
    else{
      if(mSect[ti] != tSect){
	// Try again on the other sector
	tn.SetX(cos(tPhi));
	tn.SetY(sin(tPhi));       
	tr.SetX(tRowRadius[ti]*cos(tPhi));
	tr.SetY(tRowRadius[ti]*sin(tPhi));
	// find crossing point
	tLength = helix->pathLength(tr,tn);
	tPoint = helix->at(tLength);
	if (std::isnan(tLength)){
	  cout <<"tLength loop 3rd NAN" << endl;
	 // cout <<"padrow number=  "<< ti << " not reached" << endl;
	  mSect[ti]=-1;//sector
	}
	mZ[ti] = tPoint.z();
	mSect[ti] = tSect;
	tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow, mU[ti],tPhi);
	//vector <vector <vector <TVector3 > > >tPointVec; // NEWLY ADDED >>> Vector defination to store tPoint
    	//tPointVec.push_back(tPoint);  // NEWLY ADDED >>> info storing 
	if (std::isnan(tSect)){
	  cout <<"***ERROR tSect 3"<< endl; 
	}
	if (std::isnan(tRow)){
	  cout <<"***ERROR tRow 3"<< endl;
	}
	if (std::isnan(mU[ti])){
	  cout <<"***ERROR mU[ti] 3"<< endl;
	}
	if (std::isnan(tPhi)){
	  cout <<"***ERROR tPhi 3 "<< endl;
	}  
	if(tOutOfBound || tSect!= mSect[ti] || tRow!=(ti+1)){
	  mSect[ti]=-1;
	}
      }
      }
      
      if (std::isnan(mSect[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"StHbtParticle--Fctn mSect=" << mSect[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    if (std::isnan(mU[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"StHbtParticle--Fctn mU=" << mU[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    if (std::isnan(mZ[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"StHbtParticle--Fctn mZ=" << mZ[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    // If padrow ti not reached all other beyond are not reached
    // in this case set sector to -1
    if (mSect[ti]==-1){
      for (int tj=ti; tj<45;tj++){
	mSect[tj] = -1;
	ti=45;
      }
    }
    ti++;
    if (ti<45){
     // checking helix parameters
/*cout<< " px is " << helix->momentum(-4.9845).x() << endl;
cout << "py is " << helix->momentum(-4.9845).y() << endl;
cout << "pz is " << helix->momentum(-4.9845).z() << endl;

cout << "x component of origin is " << helix->origin().x() << endl;
cout << "y component of origin is " << helix->origin().y() << endl;
cout << "z component of origin is " << helix->origin().z() << endl;

cout << "charge of the particle is " << helix->charge(-4.9845) << endl;*/
      candidates =  helix->pathLength(tRowRadius[ti]);
      tLength = (candidates.first > 0) ? candidates.first : candidates.second;}
  }
  	}
  	
 // calculate function TPClocaltransform
 int StPicoParticle::TpcLocalTransform(TVector3& aPoint, int& aSector, int& aRow,
		      float& aU, double& aPhi){
  static int tNPadAtRow[45]={
  88,96,104,112,118,126,134,142,150,158,166,174,182,
  98,100,102,104,106,106,108,110,112,112,114,116,118,120,122,122,
  124,126,128,128,130,132,134,136,138,138,140,142,144,144,144,144};
  static double tSectToPhi[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
				4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
  //static double tPhiToSect[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
	//			4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
  static double tPadWidthInner = 0.335;
  static double tPadWidthOuter = 0.67;

  static double tPi = TMath::Pi();
  // --- find sector number
  aPhi = aPoint.Phi();
  if(aPhi<0.) aPhi+=(2*tPi);
  aPhi += tPi/12.;
  if(aPhi>2*tPi) aPhi-=2*tPi;
  
  int tiPhi = (int) (aPhi/tPi*6.);
  if(aPoint.z()<0) {
    aSector = (tiPhi<3)? 3-tiPhi : 15-tiPhi;
  }
  else{
    aSector = (tiPhi<4)? 21+tiPhi : 9+tiPhi;
  }
  aPhi = tSectToPhi[aSector-1]*tPi/6.;
  //if((fabs(aPhi-aPoint.phi())>(tPi/12)){
  //cout << "Sector missmatch " << aPhi << " " << aPoint.phi() << " "
  // << aSector << endl;
  //}

  // --- calculate local coordinate
  float tR = aPoint.x()*cos(aPhi)+aPoint.y()*sin(aPhi);
  aU =      -aPoint.x()*sin(aPhi)+aPoint.y()*cos(aPhi);
 // cout <<" transformed coordinates along X axis= " << aU << endl; // NEW LINE
 // cout <<" transformed coordinates along Y axis= " << tR << endl; // NEW LINE
 
 // print global coordinates
// cout << " global x coordinates " << aPoint.x() << endl;
// cout << " global y coordinates " << aPoint.y() << endl;
  // --- find pad row 
  if(tR<57.6) {
    aRow = 0;
    return 1;
  }
  float radmax = 62.4;
  float spacing= 4.8;
  aRow=1;
  while(tR>radmax && aRow<46){
    aRow++;
    if(aRow==8){
      radmax = 96.2;
      spacing = 5.2;
    }
    else{
      if (aRow==13){
	radmax = 126.195; // lots of stuf in row 13!
	spacing = 2.0;
      }
      else{
	radmax+=spacing;
      }
    }
  }
  if(aRow>45){
    //cout << "No pad row " << tR << endl;
    return 2;
  }
  
  // --- Check if u (=aU) inbound
  double tPadWidth = aRow<14? tPadWidthInner : tPadWidthOuter;
  if(fabs(aU) > tNPadAtRow[aRow-1]*tPadWidth/2.){
    return 3;
  }

  return 0;
}	


 // calculate function to find merging pairs
 double StPicoParticle::calcMergingPar(StPicoParticle* mTrack1, StPicoParticle* mTrack2, double mMaxDuInner, double mMaxDuOuter, double mMaxDzInner, double mMaxDzOuter){
	//double mMergingParNotCalculated=0;
	double tDu, tDz;
	int tN = 0;
	double mFracOfMergedRow = 0.;
	double mWeightedAvSep =0.;
	double tDist;
	double tDistMax = 200.;
	double mClosestRowAtDCA = 0.;

	for(int ti=0 ; ti<45 ; ti++){

		if(mTrack1->mSect[ti]==mTrack2->mSect[ti] && mTrack1->mSect[ti]!=-1){
			tDu = fabs(mTrack1->mU[ti]-mTrack2->mU[ti]); 
			tDz = fabs(mTrack1->mZ[ti]-mTrack2->mZ[ti]);
			tN++;
			if(ti<13){
				mFracOfMergedRow += (tDu<mMaxDuInner && tDz<mMaxDzInner);
				tDist = ::sqrt(tDu*tDu/mMaxDuInner/mMaxDuInner+
						tDz*tDz/mMaxDzInner/mMaxDzInner);
				
			}
			else{
				mFracOfMergedRow += (tDu<mMaxDuOuter && tDz<mMaxDzOuter);
				tDist = ::sqrt(tDu*tDu/mMaxDuOuter/mMaxDuOuter+
						tDz*tDz/mMaxDzOuter/mMaxDzOuter);
				
			}
			if(tDist<tDistMax){
				mClosestRowAtDCA = ti+1;
				tDistMax = tDist;
			}
			mWeightedAvSep += tDist;
		}
	}
	if(tN>0){
		mWeightedAvSep /= tN;
		mFracOfMergedRow /= tN;
		
	}
	else{
		mClosestRowAtDCA = -1;
		mFracOfMergedRow = -1.;
		mWeightedAvSep = -1.;
		
	}
	
	return mFracOfMergedRow;
}


// Calculation of average separation
 double StPicoParticle::NominalTpcAverageSeparation(StPicoParticle* mTrack1, StPicoParticle* mTrack2) {
//	StHbtThreeVector diff;
	TVector3 diff;
	double AveSep = 0.0;
	int ipt = 0;
	if (mTrack1->tmpPosSample && mTrack2->tmpPosSample){
		while (fabs(mTrack1->tmpPosSample[ipt].x())<9999. &&
				fabs(mTrack1->tmpPosSample[ipt].y())<9999. && 
				fabs(mTrack1->tmpPosSample[ipt].z())<9999. &&
				fabs(mTrack2->tmpPosSample[ipt].x())<9999. &&
				fabs(mTrack2->tmpPosSample[ipt].y())<9999. && 
				fabs(mTrack2->tmpPosSample[ipt].z())<9999. &&
				ipt<11
		      ){
//	for (int ipt=0; ipt<11; ipt++){
//	cout << " position for track 1 " << endl;
//	mTrack1->tmpPosSample[ipt].Print();
	diff = mTrack1->tmpPosSample[ipt] - mTrack2->tmpPosSample[ipt];
	ipt++;
	AveSep += diff.Mag();
//	cout << "avg separation" << AveSep << endl;
		//	}  
	}
	AveSep = AveSep/(ipt+1.);
	return AveSep;
//	cout << "avg separation" << AveSep << endl;
	}
	else return -1;
		}


#endif

