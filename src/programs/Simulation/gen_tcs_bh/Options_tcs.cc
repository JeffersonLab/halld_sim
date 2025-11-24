#define Options_tcs_cxx
#include "Options_tcs.h"

using namespace std;

int radcor,  beamtype, mesonchoice, beampolartype, targetpol, process, model, NTotEvents, outlepton,protonorneutron,targetpoldir, beampoldir,verbose, HEP;
double Emin, E_cutoff, Emax, polbeamdeg, poltargetdeg, mt_min, mt_max, Qp2min, Qp2max, Q2max,Q2min, M_lepton ,Ebeam, Eelectron, Eproton, crossA, eta_min, eta_max, Ztarget, Atarget,targetlenght, theta_min, theta_max, M_Nucleon, thetag_max, yy_min,yy_max,philepton_min,philepton_max,Xbjmin,Xbjmax, tt, Egamma, Qp2, Q2, Phi_CMV, Theta_CMV, yy, phi_beam,Xbj, Phi_LH, mom_hms, mom_shms, theta_hms, theta_shms, thetagg_min, thetagg_max;
extern int reaction;
float param_init[30];


int Variables_TCS(){
	genSettings_t genSettings2;
	beamtype = genSettings2.beamtype;
	Emin = genSettings2.EphotonMin;
	Emax = genSettings2.EphotonMax;
	Eelectron = genSettings2.Eelectron;
	thetag_max = genSettings2.thetaphotoMax;
	NTotEvents = genSettings2.nTotEvents;
	outlepton = genSettings2.outLepton;
	targetlenght = genSettings2.targetLength;
	Atarget = genSettings2.A_target;
	Ztarget = genSettings2.Z_target;
	protonorneutron = genSettings2.protonOrNeutron;
	polbeamdeg = genSettings2.polBeamDeg;
	beampolartype = genSettings2.beamPolarType;
	targetpoldir = genSettings2.targetPolDir;
	poltargetdeg = genSettings2.polTargetDeg;
	mt_min = genSettings2.mt_Min;
	mt_max = genSettings2.mt_Max;
	Qp2min = genSettings2.Qp2Min;
	Qp2max = genSettings2.Qp2Max;
	theta_min = genSettings2.thetaCMMin;
	theta_max = genSettings2.thetaCMMax;
	Q2max = genSettings2.Q2Max;
	radcor = genSettings2.radCorrType;
	E_cutoff = genSettings2.eCut;
	HEP= genSettings2.outFormat;
	theta_hms= genSettings2.thetaHMS;
	theta_shms= genSettings2.thetaSHMS;
	mom_hms= genSettings2.pHMS;
	mom_shms= genSettings2.pSHMS;
	thetagg_min=0; thetagg_max=180;


for (int i = 0; i < 30; i++) {
    switch (i) {
        case 0: param_init[i] = (float) genSettings2.beamtype; break;
        case 1: param_init[i] = (float) genSettings2.EphotonMin; break;
        case 2: param_init[i] = (float) genSettings2.EphotonMax; break;
        case 3: param_init[i] = (float) genSettings2.Eelectron; break;
        case 4: param_init[i] = (float) genSettings2.thetaphotoMax; break;
        case 5: param_init[i] = (float) genSettings2.nTotEvents; break;
        case 6: param_init[i] = (float) genSettings2.outLepton; break;
        case 7: param_init[i] = (float) genSettings2.targetLength; break;
        case 8: param_init[i] = (float) genSettings2.A_target; break;
        case 9: param_init[i] = (float) genSettings2.Z_target; break;
        case 10: param_init[i] = (float) genSettings2.protonOrNeutron; break;
        case 11: param_init[i] = (float) genSettings2.polBeamDeg; break;
        case 12: param_init[i] = (float) genSettings2.beamPolarType; break;
        case 13: param_init[i] = (float) genSettings2.targetPolDir; break;
        case 14: param_init[i] = (float) genSettings2.polTargetDeg; break;
        case 15: param_init[i] = (float) genSettings2.mt_Min; break;
        case 16: param_init[i] = (float) genSettings2.mt_Max; break;
        case 17: param_init[i] = (float) genSettings2.Qp2Min; break;
        case 18: param_init[i] = (float) genSettings2.Qp2Max; break;
        case 19: param_init[i] = (float) genSettings2.thetaCMMin; break;
        case 20: param_init[i] = (float) genSettings2.thetaCMMax; break;
        case 21: param_init[i] = (float) genSettings2.Q2Max; break;
        case 22: param_init[i] = (float) genSettings2.radCorrType; break;
        case 23: param_init[i] = (float) genSettings2.eCut; break;
        case 24: param_init[i] = (float) genSettings2.outFormat; break;
        case 25: param_init[i] = (float) genSettings2.thetaHMS; break;
        case 26: param_init[i] = (float) genSettings2.thetaSHMS; break;
        case 27: param_init[i] = (float) genSettings2.pHMS; break;
        case 28: param_init[i] = (float) genSettings2.pSHMS; break;
        case 29: param_init[i] = (float) 0; break;  // thetagg_min
        case 30: param_init[i] = (float) 180; break;  // thetagg_max
        default: param_init[i] = 0.0f; break;
    	}
	}

	verbose=0;

	cout<<"OPTIONS::\n beamtype= "<<beamtype<<" Emin= "<<Emin<<" Emax= "<<Emax<<" Eelectron= "<<Eelectron<<endl;
	cout<<"thetag_max= "<<thetag_max<<" NTotEvents= "<<NTotEvents<<" outlepton= "<<outlepton<<endl;
	cout<<"targetlenght= "<<targetlenght<<" Atarget= "<<Atarget<<" Ztarget= "<<Ztarget<<" protonorneutron "<<protonorneutron<<endl;
	cout<<"polbeamdeg= "<<polbeamdeg<<" targetpoldir= "<<targetpoldir<<" poltargetdeg= "<<poltargetdeg<<endl;
	cout<<"mt_min= "<<mt_min<<" mt_max= "<<mt_max<<" Qp2min= "<<Qp2min<<" Qp2max= "<<Qp2max<< " Q2max= "<<Q2max<<endl;
	cout<<"theta_min= "<<theta_min<<" theta_max= "<<theta_max<<endl;
	cout<<"data output= "<<HEP<<endl;
	if (HEP>=3) cout<<"hms/shms limits "<<theta_hms<<" "<<mom_hms<<" "<<theta_shms<<" "<<mom_shms<<endl;
	return 1;
}


int Variables_PSEEPHOTO_FIX(){
	genSettings_t genSettings2;
	beamtype = genSettings2.beamtype;
	Emin = genSettings2.EphotonMin;
	Emax = genSettings2.EphotonMax;
	Eelectron = genSettings2.Eelectron;
	thetag_max = genSettings2.thetaphotoMax;
	NTotEvents = genSettings2.nTotEvents;
	outlepton = genSettings2.outLepton;
	targetlenght = genSettings2.targetLength;
	Atarget = genSettings2.A_target;
	Ztarget = genSettings2.Z_target;
	mt_min = genSettings2.mt_Min;
	mt_max = genSettings2.mt_Max;
	Qp2min = genSettings2.Qp2Min;
	Qp2max = genSettings2.Qp2Max;
	theta_min = genSettings2.thetaCMMin;
	theta_max = genSettings2.thetaCMMax;
	Q2max = genSettings2.Q2Max;
	radcor = genSettings2.radCorrType;
	E_cutoff = genSettings2.eCut;
	HEP= genSettings2.outFormat;
	theta_hms= genSettings2.thetaHMS;
	theta_shms= genSettings2.thetaSHMS;
	mom_hms= genSettings2.pHMS;
	mom_shms= genSettings2.pSHMS;
	thetagg_min=0; thetagg_max=180;


for (int i = 0; i < 25; i++) {
    switch (i) {
        case 0: param_init[i] = (float) genSettings2.beamtype; break;
        case 1: param_init[i] = (float) genSettings2.EphotonMin; break;
        case 2: param_init[i] = (float) genSettings2.EphotonMax; break;
        case 3: param_init[i] = (float) genSettings2.Eelectron; break;
        case 4: param_init[i] = (float) genSettings2.thetaphotoMax; break;
        case 5: param_init[i] = (float) genSettings2.nTotEvents; break;
        case 6: param_init[i] = (float) genSettings2.outLepton; break;
        case 7: param_init[i] = (float) genSettings2.targetLength; break;
        case 8: param_init[i] = (float) genSettings2.A_target; break;
        case 9: param_init[i] = (float) genSettings2.Z_target; break;
        case 10: param_init[i] = (float) genSettings2.mt_Min; break;
        case 11: param_init[i] = (float) genSettings2.mt_Max; break;
        case 12: param_init[i] = (float) genSettings2.Qp2Min; break;
        case 13: param_init[i] = (float) genSettings2.Qp2Max; break;
        case 14: param_init[i] = (float) genSettings2.thetaCMMin; break;
        case 15: param_init[i] = (float) genSettings2.thetaCMMax; break;
        case 16: param_init[i] = (float) genSettings2.Q2Max; break;
        case 17: param_init[i] = (float) genSettings2.radCorrType; break;
        case 18: param_init[i] = (float) genSettings2.eCut; break;
        case 19: param_init[i] = (float) genSettings2.outFormat; break;
        case 20: param_init[i] = (float) genSettings2.thetaHMS; break;
        case 21: param_init[i] = (float) genSettings2.thetaSHMS; break;
        case 22: param_init[i] = (float) genSettings2.pHMS; break;
        case 23: param_init[i] = (float) genSettings2.pSHMS; break;
        case 24: param_init[i] = (float) 0; break;  // thetagg_min
        case 25: param_init[i] = (float) 180; break;  // thetagg_max
        default: param_init[i] = 0.0f; break;
		}
	}

	verbose=0;

	cout<<"OPTIONS::\n beamtype= "<<beamtype<<" Emin= "<<Emin<<" Emax= "<<Emax<<" Eelectron= "<<Eelectron<<endl;
	cout<<"thetag_max= "<<thetag_max<<" NTotEvents= "<<NTotEvents<<" outlepton= "<<outlepton<<endl;
	cout<<"targetlenght= "<<targetlenght<<" Atarget= "<<Atarget<<" Ztarget= "<<Ztarget<<" protonorneutron "<<protonorneutron<<endl;
	cout<<"polbeamdeg= "<<polbeamdeg<<" targetpoldir= "<<targetpoldir<<" poltargetdeg= "<<poltargetdeg<<endl;
	cout<<"mt_min= "<<mt_min<<" mt_max= "<<mt_max<<" Qp2min= "<<Qp2min<<" Qp2max= "<<Qp2max<< " Q2max= "<<Q2max<<endl;
	cout<<"theta_min= "<<theta_min<<" theta_max= "<<theta_max<<endl;
	cout<<"output format: "<<HEP<<endl;
	if (HEP>=3) cout<<"hms/shms limits "<<theta_hms<<" "<<mom_hms<<" "<<theta_shms<<" "<<mom_shms<<endl;
	return 1;
}
