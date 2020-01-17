/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
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
#include "MyReadConfig.h"

int main(int argc, char * argv[]) {
  
  Int_t NbOfEvt = 0;
  Char_t * genconfigfile = 0;
  Int_t run_nb_min = 0;
  Int_t run_nb_max = 0;
  Char_t * workflow = 0;
  Char_t * out_dir = 0;
  if (argc == 7) {
    NbOfEvt = atof(argv[1]);
    genconfigfile = argv[2];
    run_nb_min = atof(argv[3]);
    run_nb_max = atof(argv[4]);
    workflow = argv[5];
    out_dir = argv[6];
  } else {
    std::cout << "run_whizard" << std::endl;
    std::cout << "===========" << std::endl;
    std::cout << "This program calculates the number of events to be launched for different fixed incident photon-beam energy" << std::endl;
    std::cout << "                                                                 , launches WHIZAR, and produced a LHE file"<< std::endl;
    std::cout << std::endl;
    std::cout << "Usage:                                                                                                     " << std::endl;
    std::cout << "Method : ./run_whizard NbOfEvt generator_configfile RunNbMin RunNbMax workflow_name output_directory_path  " << std::endl;
  }
  gRandom->SetSeed(0);
  
  // Get generator config file
  MyReadConfig * ReadWOFile = new MyReadConfig();
  ReadWOFile->ReadConfigFile(genconfigfile);
  TString m_run_wo = ReadWOFile->GetConfigName("run_wo"); 
  TString m_process = ReadWOFile->GetConfigName("process"); 
  TString m_lhe_dir = ReadWOFile->GetConfigName("lhe_dir"); 
  TString m_rootfile = ReadWOFile->GetConfigName("ROOTFluxFile"); 
  TString m_histoname = ReadWOFile->GetConfigName("ROOTFluxName");
  Double_t * beamLowE = ReadWOFile->GetConfig1Par("PhotonBeamLowEnergy");
  Double_t * beamHighE = ReadWOFile->GetConfig1Par("PhotonBeamHighEnergy");
  cout << " m_run_wo " << m_run_wo << " beamLowE " << beamLowE[0] << " GeV beamHighE " << beamHighE[0] << endl;
  
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
  
  TFile * diagOut = new TFile( TString::Format("gen_whizard_%s_run_nb_min_%d_max_%d.root", m_process.Data(), run_nb_min, run_nb_max), "recreate" );
  TH1F * h_egam = new TH1F("egam", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  
  if (m_run_wo == "true") {
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
      egam *= 1e3;
      int nbofevt =  h_egam->GetBinContent(i + 1);
      if (nbofevt > 0) {
	int seed_nb = gRandom->Uniform(0, 1e9);
	system(TString::Format("./whizard.sh %s %d %d %d %d %d %s %s", m_process.Data(), nbofevt, (int) egam, run_nb_min, run_nb_max, seed_nb, workflow, out_dir));
	npts ++;
      }
    }
  }
  h_egam->Write();
  diagOut->Close();
}
