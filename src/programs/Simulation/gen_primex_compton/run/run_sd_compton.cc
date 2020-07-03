/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       HallD Group                                     * 
*                                                                         *                                                                                                                               
* Author: The HallD Group                                                 *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include "Riostream.h"
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include "TRandom3.h"
#include "UTILITIES/MyReadConfig.h"

int main(int argc, char * argv[]) {
  
  Int_t NbOfEvt = 0;
  Char_t * genconfigfile = 0;
  Int_t run_nb_min = 0;
  Int_t run_nb_max = 0;
  Char_t * workflow = 0;
  Char_t * out_dir = 0;
  Char_t * shell = 0;
  int iswf = 1;
  if (argc == 7) {
    NbOfEvt = atof(argv[1]);
    genconfigfile = argv[2];
    run_nb_min = atof(argv[3]);
    run_nb_max = atof(argv[4]);
    iswf = 0;
    out_dir = argv[5];
    shell = argv[6];
  } else if (argc == 8) {
    NbOfEvt = atof(argv[1]);
    genconfigfile = argv[2];
    run_nb_min = atof(argv[3]);
    run_nb_max = atof(argv[4]);
    workflow = argv[5];
    out_dir = argv[6];
    shell = argv[7];
  } else {
    std::cout << "run_whizard" << std::endl;
    std::cout << "===========" << std::endl;
    std::cout << "This program calculates the number of events to be launched for different fixed incident photon-beam energy, launches  " << std::endl;
    std::cout << "                                                                                 sd_comptom, and produced a root file  " << std::endl;
    std::cout << std::endl;
    std::cout << "Usage:                                                                                                                 " << std::endl;
    std::cout << "Method 1 :run_sd_compton NbOfEvt generator_configfile RunNbMin RunNbMax output_directory_path your_shell               " << std::endl;
    std::cout << "Method 2 :run_sd_compton NbOfEvt generator_configfile RunNbMin RunNbMax workflow_name output_directory_path your_shell " << std::endl;
  }
  gRandom->SetSeed(0);
  
  // Get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(genconfigfile);
  TString m_rootfile = ReadFile->GetConfigName("ROOTFluxFile"); 
  TString m_histoname = ReadFile->GetConfigName("ROOTFluxName");
  Double_t * beamLowE = ReadFile->GetConfig1Par("PhotonBeamLowEnergy");
  Double_t * beamHighE = ReadFile->GetConfig1Par("PhotonBeamHighEnergy");
  cout << " beamLowE " << beamLowE[0] << " GeV beamHighE " << beamHighE[0] << endl;
  
  if (m_histoname == "FLUXNAME") 
    m_rootfile = "";

  // Assume a beam energy spectrum of 1/E(gamma)
  TF1 ebeam_spectrum("beam_spectrum","1/x",beamLowE[0], beamHighE[0]);
  
  // Get beam properties from configuration file
  TFile * ifile = new TFile(m_rootfile);
  TH1F * cobrem_vs_E = (TH1F *) ifile->Get(m_histoname);
  for (int i = 0; i <  cobrem_vs_E->GetNbinsX(); i ++) {
    double egam = cobrem_vs_E->GetBinCenter(i + 1);
    if (egam < beamLowE[0] || egam > beamHighE[0]) 
      cobrem_vs_E->SetBinContent(i + 1, 0.);
  }
  
  TFile * diagOut = new TFile( TString::Format("gen_sd_Compton_run_nb_min_%d_max_%d.root", run_nb_min, run_nb_max), "recreate" );
  TH1F * h_egam = new TH1F("egam", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  
  for (int i = 0; i < NbOfEvt; ++i) {
    if (i%10000 == 1)
      cout << "event " << i <<endl;
    
    // get beam energy
    double ebeam = 0;
    if (m_rootfile == "" || cobrem_vs_E == 0) 
      ebeam = ebeam_spectrum.GetRandom();
    else if (m_rootfile != 0)
      ebeam = cobrem_vs_E->GetRandom();
    
    h_egam->Fill(ebeam);
  }
  //cout << " bin nb " << h_egam->GetNbinsX() << endl;
  int npts = 0;
  
  for (int i = 0; i < h_egam->GetNbinsX(); i ++) {
    double egam = h_egam->GetBinCenter(i + 1);
    int nbofevt =  h_egam->GetBinContent(i + 1);
    if (nbofevt > 0) {
      //egam *= 1e3;
      if (strcmp(shell,"tcsh") == 0 && iswf == 0)
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d", 
			       //out_dir, (int) egam, nbofevt, run_nb_min, run_nb_max));
			       out_dir, egam, nbofevt, run_nb_min, run_nb_max));
      else if (strcmp(shell,"tcsh") == 0 && iswf == 1)
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.csh %s %f %d %d %d %s", 
			       //out_dir, (int) egam, nbofevt, run_nb_min, run_nb_max, workflow));
			       out_dir, egam, nbofevt, run_nb_min, run_nb_max, workflow));
      else if (strcmp(shell,"bash") == 0 && iswf == 0)
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d", 
			       //out_dir, (int) egam, nbofevt, run_nb_min, run_nb_max));
			       out_dir, egam, nbofevt, run_nb_min, run_nb_max));
      else if (strcmp(shell,"bash") == 0 && iswf == 1)
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.sh %s %f %d %d %d %s", 
			       //out_dir, (int) egam, nbofevt, run_nb_min, run_nb_max, workflow));
			       out_dir, egam, nbofevt, run_nb_min, run_nb_max, workflow));
      npts ++;
    }
  }
  
  h_egam->Write();
  diagOut->Close();
}
