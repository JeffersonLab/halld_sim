#ifndef TreatOptions_h
#define TreatOptions_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Options_tcs.h"
#include "Utils.h"
#include <unistd.h>
#include "TableProvider.h"
#include "FormFactors.h"
#include "RadProcess.h"
#include "FillTables.h"
#include "Constants.h"
#include "RDISparam.h"
#include "PDFparam.h"
#include "genSettings_t.h"


    // Binning for TCS
    extern int NEB, NT, NQp2, NTh, NPhi, NPhis, NPsis; // tablev11 costh<1
    extern float CosThetaMaxTCS;// 0.9995;
    extern double EMINI,EMAXI,TMINI,TMAXI,QP2MINI,QP2MAXI; // table v11
    extern double PHIMINI, PHIMAXI, THMINI, THMAXI; 
    extern double PHISMINI, PHISMAXI, PSISMINI, PSISMAXI;

    // table cuts (TCS)
    extern int thmE, thmQ, thmT;
    extern float thmEmin,thmEmax,thmQmin,thmQmax,thmTmin,thmTmax;  
   
    void settcsbins(int Qp2_LimitType3);
    int VariablesValues_TCS();
    int VariablesValues_PSTCSFIX();
    int RandomGen_TCS();

    double BH_anal_TCS_exact(double eg, double q_out_th_lab, double pair_e_th_cm, double pair_e_phi_cm);

    //int FillTable(int reaction, TableProvider& table_provider);

     double L_BH_TCS(double eg, double q_out_th_lab, double pair_e_th_cm, double pair_e_phi_cm); 

// ****** High Q'2 TCS ******
	// Binning for TCS
/*    const int NEB=13, NT=33, NQp2=18, NTh=20, NPhi=20, NPhis=20, NPsis=20; // tablev11 costh<1
    const float CosThetaMaxTCS = 1;// 0.9995; 
    const double EMINI=5.0,EMAXI=11.5,TMINI=0.04,TMAXI=2.02,QP2MINI=3.8,QP2MAXI=9.2; // table v11
    const double PHIMINI=3., PHIMAXI=363., THMINI=30., THMAXI=150., PHISMINI = 0, PHISMAXI=360, PSISMINI = 0, PSISMAXI = 360;

    // table cuts (TCS)
    const int thmE=14, thmQ=20, thmT=40;
    const float thmEmin=5,thmEmax=12,thmQmin=3.8,thmQmax=9.3,thmTmin=0.02,thmTmax=2.02;
*/

// ****** Low Q'2 TCS ******
	// Binning for TCS
/*   const int NEB=13, NT=30, NQp2=15, NTh=25, NPhi=20, NPhis=20, NPsis=20; // tablev11 costh<1
    const float CosThetaMaxTCS = 1;// 0.9995;
    const double EMINI=5.0,EMAXI=11.5,TMINI=0.04,TMAXI=1.54,QP2MINI=0.8,QP2MAXI=5.3; // table v11
    const double PHIMINI=3., PHIMAXI=363., THMINI=30., THMAXI=155., PHISMINI = 0, PHISMAXI=360, PSISMINI = 0, PSISMAXI = 360;

    // table cuts (TCS)
    const int thmE=14, thmQ=30, thmT=40;
    const float thmEmin=5,thmEmax=12,thmQmin=0.8,thmQmax=9.3,thmTmin=0.02,thmTmax=2.02; 
*/

#endif


