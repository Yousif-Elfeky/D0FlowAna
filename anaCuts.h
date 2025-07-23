const int ncent = 4;
const char nameCent[ncent][250] = {"0-10%","10-40%","40-80%","0-80%"};
const char nameCent1[ncent][250] = {"0_10","10_40","40_80","0_80"};
const float centLw[ncent] = {7,4,0,0}; // >=
const float centUp[ncent] = {9,7,4,9}; // <

const float D0BR = 0.0388;
//const int effID[ncent] = {2, 3, 4, 1};
//const float Nev[ncent] = {9.82e7, 2.75e8, 3.71e8, 7.45e8};


const int npt = 8;
const double nptbin[npt+1] = {0.,.7,1.1,1.6,2.2,3.,4.,5.,8.};

//const int npt = 11;
//const double nptbin[npt+1] = { 0., 0.5, 1., 1.5, 2., 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0 };

const int npt1 = 5;
const double nptbin1[npt+1] = {0., 1.4, 2.0, 2.5, 3.0, 5.};

//float fitRange_lw[npt] = {1.725, 1.725,1.725, 1.725, 1.725, 1.725};
//float fitRange_up[npt] = {2.025, 2.025,2.025, 2.025, 2.025, 2.025};
float fitRange_lw[npt] = {1.69, 1.69,1.69, 1.69,1.69, 1.69, 1.69, 1.69};
float fitRange_up[npt] = {2.04, 2.04,2.04, 2.04,2.04, 2.04, 2.04, 2.04};

float norm_tune[4] = {0.101, 0.10164, 0.1025, 0.1034};

// mix event scale range (to Like sign)

//const float mR_lw = 1.725;
//const float mR_up = 2.025;
const float mR_lw = 1.69;
const float mR_up = 2.0;


const float inc = 0.00002;
const float normScan_lw = -0.0003;
const float normScan_up = 0.0003;
