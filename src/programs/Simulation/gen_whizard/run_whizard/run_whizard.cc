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

  Int_t NbOfEvt = atof(argv[1]);
  cout <<"NbOfEvt " << NbOfEvt << endl;
  //Char_t * beamconfigfile = argv[2];
  //cout << "beamconfigfile " << beamconfigfile << endl;
  Char_t * genconfigfile = argv[2];
  cout << "genconfigfile " << genconfigfile << endl;
  Int_t run_nb = atof(argv[3]);
  gRandom->SetSeed(0);
  int seed_nb = gRandom->Uniform(0, 1e9);
  Char_t * workflow = argv[4];
  Char_t * out_dir = argv[5];
  
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

  TFile * diagOut = new TFile( TString::Format("gen_whizard_%s.root", m_process.Data()), "recreate" );
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
	//cout << m_process << " nbofevt " <<  nbofevt << " egam " << (int) egam << " seed " << seed_nb << endl;
        system(TString::Format("./whizard.sh %s %d %d %d %d %s %s", m_process.Data(), nbofevt, (int) egam, run_nb, seed_nb, workflow, out_dir));
	//cout << "egam " << egam << endl;
	npts ++;
      }
    }
    //cout << "Number of different photon energies thrown " << npts << endl;
  }
  h_egam->Write();
  diagOut->Close();
}
