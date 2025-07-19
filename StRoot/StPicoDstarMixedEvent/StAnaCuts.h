#ifndef StAnaCuts_H
#define StAnaCuts_H
#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
    std::array<unsigned int, 4> const triggers = {600001, 600011, 600021, 600031};
     //event level
    float const Vr = 2.0; //cm
    float const Vz = 30.0;
	float const VzMax = 25.0;
	float const VzMin = -35.0;
    float const VzDiff = 3.;
    float const Verror = 1.0e-5;
    //track level
    float const minPtCut1 = 0.15;
    float const maxPtCut1 = 2.;
    float const EtaCut = 1.;
    int const NHitsFit1 = 15;
	
    float const minPtCut2 = 0.6;
    int const NHitsFit2 = 20;
    
	float const PtsRMin = 0.52;
    float const PtsRMax = 1.2;
    float const Dca1 = 3.;
    float const Dca2 = 2.;
	float const mNHitsDedx = 10;
    float const PIDPpCut = 1.6;
    //D0
    float const D0yCut = 1.;

    
}
#endif
