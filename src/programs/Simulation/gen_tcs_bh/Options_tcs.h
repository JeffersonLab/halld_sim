#ifndef Options_tcs_h
#define Options_tcs_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "genSettings_t.h"

int Variables_TCS();
int Variables_PSEEPHOTO_FIX();

//double VFluxFactor(double Q2_max, double E, double nu);
extern int mesonchoice, radcor, beamtype, beampolartype, targetpol, process, model, NTotEvents, outlepton,protonorneutron, targetpoldir, beampoldir,verbose, HEP ;
extern double E_cutoff, Ebeam, Emin, Emax, polbeamdeg, poltargetdeg, mt_min, mt_max, Qp2min, Qp2max, Q2max,Q2min, M_lepton,  Eelectron, Xbjmin,Xbjmax, thetag_max, Eproton, crossA, eta_min, eta_max, Ztarget, Atarget,targetlenght, theta_min, theta_max, M_Nucleon, yy_min,yy_max,philepton_min,philepton_max, tt, Egamma, Qp2, Q2, Phi_CMV, Theta_CMV, yy, phi_beam,Xbj, Phi_LH, theta_hms,theta_shms, mom_hms, mom_shms, thetagg_min, thetagg_max;
extern float param_init[30];

#endif
