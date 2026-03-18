#define TreeInit_cxx
#include <iostream>
#include <fstream>
#include "TreeInit.h"
#include <limits>
#include <unistd.h>
using namespace std;

void TreeInit (TTree *Dump_Tree){ //, int process, int pol ){

/*        Dump_Tree->Branch("indexrun",&indexrun, "indexrun/I");
        Dump_Tree->Branch("param_initfile", &param_initfile,"param_initfile[30]/F");
        Dump_Tree->Branch("phase_space", &phase_space, "phase_space/F");
        Dump_Tree->Branch("ALV_minus_lab",&ALV_minus_lab,"ALV_minus_lab[4]/D");
        Dump_Tree->Branch("ALV_plus_lab",&ALV_plus_lab,"ALV_plus_lab[4]/D");
        Dump_Tree->Branch("ALV_gamma_in",&ALV_gamma_in,"ALV_gamma_in[4]/D");
        Dump_Tree->Branch("ALV_gamma_out_lab",&ALV_gamma_out_lab,"ALV_gamma_out_lab[4]/D");
        Dump_Tree->Branch("ALV_el_in",&ALV_el_in,"ALV_el_in[4]/D");
        Dump_Tree->Branch("ALV_Recoil_lab",&ALV_Recoil_lab,"ALV_Recoil_lab[4]/D");
        Dump_Tree->Branch("ALV_el_out",&ALV_el_out,"ALV_el_out[4]/D");
        Dump_Tree->Branch("ALV_minus_CMV",&ALV_minus_CMV,"ALV_minus_CMV[4]/D");
        Dump_Tree->Branch("ALV_plus_CMV",&ALV_plus_CMV,"ALV_plus_CMV[4]/D");
        Dump_Tree->Branch("ALV_minus_CMeP",&ALV_minus_CMeP,"ALV_minus_CMeP[4]/D");
        Dump_Tree->Branch("ALV_plus_CMeP",&ALV_plus_CMeP,"ALV_plus_CMeP[4]/D");
        Dump_Tree->Branch("ALV_gamma_CMeP",&ALV_gamma_CMeP, "ALV_gamma_CMeP[4]/D");
        Dump_Tree->Branch("ALV_Recoil_CMeP",&ALV_Recoil_CMeP,"ALV_Recoil_CMeP[4]/D");
        Dump_Tree->Branch("ALV_Virtual_CMeP",&ALV_Virtual_CMeP, "ALV_Virtual_CMeP[4]/D");
        Dump_Tree->Branch("ALV_Virtual",&ALV_Virtual, "ALV_Virtual[4]/D");
        Dump_Tree->Branch("ALV_Target_CMeP",&ALV_Target_CMeP,"ALV_Target_CMeP[4]/D");
        Dump_Tree->Branch("Q2",&Q2,"Q2/D");
        Dump_Tree->Branch("theta_beam",&theta_beam ,"theta_beam/D");
        Dump_Tree->Branch("phi_beam",&phi_beam ,"phi_beam/D");
  */  
    Dump_Tree->Branch("yy",&yy ,"yy/D");
        Dump_Tree->Branch("Xbj",&Xbj ,"Xbj/D");
        Dump_Tree->Branch("Qp2",&Qp2,"Qp2/D");
        Dump_Tree->Branch("tt",&tt,"tt/D");
    /*    Dump_Tree->Branch("ttmin",&ttmin,"ttmin/D");
        Dump_Tree->Branch("Phi_CMV",&Phi_CMV,"Phi_CMV/D");
        Dump_Tree->Branch("Phi_LH",&Phi_LH,"Phi_LH/D");
        Dump_Tree->Branch("Theta_CMV",&Theta_CMV,"Theta_CMV/D");
        Dump_Tree->Branch("EventNumber", &EventNumber , "EventNumber/L");
        Dump_Tree->Branch("TrueEventNumber", &TrueEventNumber , "TrueEventNumber/L");
*/
}




