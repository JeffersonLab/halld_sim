/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstdlib>

#include "particleType.h"

#include "UTILITIES/BeamProperties.h"

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"

#include "AMPTOOLS_AMPS/Compton.h"

#include "AMPTOOLS_MCGEN/GammaPToXP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"
#include "TSystem.h"
#include <TSystemDirectory.h>
#include "Riostream.h"
#include "TGraph.h"
#include "TChain.h"

#include "HddmOut.h"
#include "MyReadConfig.h"
#include "HDDM/hddm_s.hpp"

using std::complex;
using namespace std;

#define GAMMA_TYPE 1
#define ELECTRON_TYPE 3

int main( int argc, char* argv[] ){
  
  string  beamconfigfile("");
  TString genconfigfile("");// = "gen_config.dat";
  TString whizard("");
  string  outname("");
  string  hddmname("");
  
  ifstream in_coherent;
  ifstream in_incoherent;
  
  double beamLowE   = 3.0;
  double beamHighE  = 12.0;
  const double M_electron = 0.000511;
  const double M_gamma = 0.0;
  TLorentzVector Target_4Vec(0, 0, 0, M_electron);
  
  int runNum = 9001;
  int seed = 0;
  
  int nEvents = 10000;
  double Na = 6.002214199 * 1e23; // mol^-1
  double alpha = 1.0 / 137.036;
  double bohrRadius = (1.0 / alpha) * (1.0 / M_electron);
  
  //parse command line:
  for (int i = 1; i < argc; i++) {
    
    string arg(argv[i]);
    
    
    if (arg == "-sa") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  whizard = argv[++i]; 
    }
    if (arg == "-c") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  beamconfigfile = argv[++i]; 
    }
    if (arg == "-e") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  genconfigfile = argv[++i]; 
    }
    if (arg == "-o") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outname = argv[++i]; 
    }
    if (arg == "-hd") {
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  hddmname = argv[++i]; 
    }
    if (arg == "-n"){   
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  nEvents = atoi( argv[++i] ); 
    }
    if (arg == "-a") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  beamLowE = atof( argv[++i] ); 
    }
    if (arg == "-b") {  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  beamHighE = atof( argv[++i] ); 
    }
    if (arg == "-r") {
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  runNum = atoi( argv[++i] ); 
    }
    if (arg == "-s") {
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  seed = atoi( argv[++i] ); 
    }
    if (arg == "-h") {
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -c  <file>\t Beam config file" << endl;
      cout << "\t -e  <file>\t Generator config file" << endl;
      cout << "\t -hd <name>\t HDDM file output name [optional]" << endl;
      cout << "\t -n  <value>\t Number of events to generate [optional]" << endl;
      cout << "\t -a  <value>\t Minimum photon energy to simulate events [optional]" << endl;
      cout << "\t -b  <value>\t Maximum photon energy to simulate events [optional]" << endl;
      cout << "\t -r  <value>\t Run number assigned to generated events [optional]" << endl;
      cout << "\t -s  <value>\t Random number seed initialization [optional]" << endl;
      exit(1);
    }
  }
  
  if (outname.size() == 0 && hddmname == 0) {
    cout << "No output specificed:  run gen_compton_simple -h for help" << endl;
    exit(1);
  }
  
  if (genconfigfile == "") {
    cout << "No generator configuration file: run gen_primex_eta_he4 -h for help " << endl;
    exit(1);
  }
  // random number initialization (set to 0 by default)
  gRandom->SetSeed(seed);
  
  int seed_nb = gRandom->Uniform(0, 1e9);

  // initialize HDDM output
  HddmOut *hddmWriter = nullptr;
  if (hddmname != "")
    hddmWriter = new HddmOut(hddmname.c_str());
  
  // Assume a beam energy spectrum of 1/E(gamma)
  TF1 ebeam_spectrum("beam_spectrum","1/x",beamLowE,beamHighE);
  
  // Get beam properties from configuration file
  TH1D * cobrem_vs_E = 0;
  if (beamconfigfile != "") {
    BeamProperties beamProp( beamconfigfile );
    cobrem_vs_E = (TH1D*)beamProp.GetFlux();
  }
  
  // Get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(genconfigfile);
  TString m_workflow = ReadFile->GetConfigName("worflow"); 
  TString m_run_wo = ReadFile->GetConfigName("run_wo"); 
  TString m_process = ReadFile->GetConfigName("process"); 
  TString m_lhe_dir = ReadFile->GetConfigName("lhe_dir"); 
  TString m_lhe_file = ReadFile->GetConfigName("lhe_file"); 
  TString m_out_dir = ReadFile->GetConfigName("out_dir"); 
  Double_t * m_target = ReadFile->GetConfig4Par("target");  
  TString m_XS_pair = ReadFile->GetConfigName("XS_pair"); 
  TString m_XS_trip = ReadFile->GetConfigName("XS_trip"); 
  double Z = m_target[0];
  double A = m_target[1];
  double rho = m_target[2];
  double Ltarget = m_target[3];
  
  double Luminosity = Na / A * rho * Ltarget * 1e-27; // cm^2 to mb^-1

  TFile * diagOut = new TFile( TString::Format("gen_whizard_%s.root", m_process.Data()), "recreate" );
  TH1F * h_egam1 = new TH1F("egam1", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_egam2 = new TH1F("egam2", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_Tkin_gam = new TH1F("Tkin_gam", ";T_{#gamma}^{kin} [GeV];Count/10MeV", 1200, 0.0, 12.0);
  TH1F * h_Tkin_rec = new TH1F("Tkin_rec", ";T_{e^{-}-recoil}^{kin} [GeV];Count/10MeV", 1200, 0.0, 12.0);
  TH2F * h_theta_vs_Tkin_gam = new TH2F("theta_vs_Tkin_gam", ";T_{#gamma}^{kin} [GeV];log_{10}(#theta) [^{o}];Count/10MeV", 1200, 0.0, 12.0, 1200, -10, 2.25);
  TH2F * h_theta_vs_Tkin_rec = new TH2F("theta_vs_Tkin_rec", ";T_{e^{-}-recoil}^{kin} [GeV];log_{10}(#theta) [^{o}];Count/10MeV", 1200, 0.0, 12.0, 1200, -10, 2.25);
  TH1F * h_lgam1 = new TH1F("lgam1", ";E_{#gamma} [GeV]; Luminosity MeV^{1} #cdot mb^{-1}", 12000, 0.0, 12.0);
  TH1F * h_lgam2 = new TH1F("lgam2", ";E_{#gamma} [GeV]; Luminosity MeV^{1} #cdot mb^{-1}", 12000, 0.0, 12.0);
  TH1F * h_lgam3 = new TH1F("lgam3", ";E_{#gamma} [GeV]; Luminosity MeV^{1} #cdot mb^{-1}", 12000, 0.0, 12.0);
  
  if (m_run_wo == "true") { //If "true" run WO
    
    for (int i = 0; i < nEvents; ++i) { //Generate photon beam spectrum
      if (i%10000 == 1)
	cout << "event " << i <<endl;
      
      // get beam energy
      double ebeam = 0;
      if (beamconfigfile == "" || cobrem_vs_E == 0) 
	ebeam = ebeam_spectrum.GetRandom();
      else if (beamconfigfile != "")
	ebeam = cobrem_vs_E->GetRandom();
      
      h_egam1->Fill(ebeam);
    }
    
    for (int i = 0; i < h_egam1->GetNbinsX(); i ++) { //Generate LHE file
      double egam = h_egam1->GetBinCenter(i + 1);
      egam *= 1e3;
      int nbofevt =  h_egam1->GetBinContent(i + 1);
      if (nbofevt > 0) {
	system(TString::Format("./whizard.sh %s %d %d %d %d %s %s", m_process.Data(), nbofevt, (int) egam, runNum, seed_nb, m_workflow.Data(), m_out_dir.Data()));
      }
    }
  }
  
  if (m_run_wo == "false" || m_lhe_dir != "" || m_lhe_file != "") { //Read and loop over a single or all root file
    
    TSystemDirectory dir(m_lhe_dir.Data(), m_lhe_dir.Data());
    TList * files = dir.GetListOfFiles(); 
    
    TGraph * grXS_pair = new TGraph(m_XS_pair);
    TGraph * grXS_trip = new TGraph(m_XS_trip);
    
    int npart = 0;
    int pdg[100];
    double xs = 0, er_xs = 0, weight = 0;
    double px[100], py[100], pz[100], e[100], m[100];
    int status[100], first_daughter[100], last_daughter[100];
    for (int i = 0; i < 100; i ++) {
      pdg[i] = 0;
      status[i] = 0;
      first_daughter[i] = 0;
      last_daughter[i] = 0;
      px[i] = 0;
      py[i] = 0;
      pz[i] = 0;
      e[i] = 0;
      m[i] = 0;
    }
    
    TChain * m_tree = new TChain("lhe");
    
    if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next()) ) {
	fname = file->GetName();
	if (!file->IsDirectory() && fname.Contains(".root") && fname.Contains(m_process)) {
	  TString RootFileName = m_lhe_dir + fname;
	  m_tree->Add(RootFileName);
	}
      }
    }
    
    if (m_lhe_file != "") 
      m_tree->Add(m_lhe_file);
    
    Long64_t NbOfEvent = m_tree->GetEntries();
    
    m_tree->SetBranchAddress("npart",&npart);
    m_tree->SetBranchAddress("weight",&weight);
    m_tree->SetBranchAddress("xs",&xs);
    m_tree->SetBranchAddress("er_xs",&er_xs);
    m_tree->SetBranchAddress("pdg", pdg);
    m_tree->SetBranchAddress("status", status);
    m_tree->SetBranchAddress("first_daughter", first_daughter);
    m_tree->SetBranchAddress("last_daughter", last_daughter);
    m_tree->SetBranchAddress("px", px);
    m_tree->SetBranchAddress("py", py);
    m_tree->SetBranchAddress("pz", pz);
    m_tree->SetBranchAddress("e", e);
    m_tree->SetBranchAddress("m", m);  
    
    for (int counter = 0; counter < NbOfEvent; counter ++) {
      
      m_tree->GetEntry(counter);
      
      TLorentzVector gamma_4Vec(0, 0, 0, 0);
      TLorentzVector photon_4Vec[10];
      TLorentzVector electron_4Vec[10];
      int npart_photon = 0;
      int npart_electron = 0;
      double Emax = 0;
      int nid_recoil = 0;
      TLorentzVector e_recoil_4Vec(0,0,0,0);
      for (int i = 0; i < npart; i ++) {
	photon_4Vec[i] = TLorentzVector(0,0,0,0);
	electron_4Vec[i] = TLorentzVector(0,0,0,0);
	if (status[i] == -1 && pdg[i] == 22) gamma_4Vec = TLorentzVector(0, 0, e[i], e[i]);
	if (status[i] == 1 && pdg[i] == 11) {
	  electron_4Vec[npart_electron] = TLorentzVector(px[i], py[i], pz[i], e[i]);
	  if (e[i] > Emax) {
	    Emax = e[i];
	    e_recoil_4Vec = TLorentzVector(px[i], py[i], pz[i], e[i]);
	    nid_recoil = i;
	  }
	  npart_electron ++;
	}
	if (status[i] == 1 && pdg[i] == 22) {
	  photon_4Vec[npart_photon] = TLorentzVector(px[i], py[i], pz[i], e[i]); 
	  npart_photon ++;
	}
      }
      
      TLorentzVector moTransfer = e_recoil_4Vec - Target_4Vec;
      double e_gamma = gamma_4Vec.E();
      h_egam2->Fill(e_gamma);
      h_lgam1->Fill(e_gamma, Luminosity);
      //Calculation of screening and radiative corrections
      //Screening factor
      double SF = 1.0 / pow(1 + pow(bohrRadius * moTransfer.P() / 2.0, 2), 2);
      double ScreeningFactor = 1.0 - pow(SF, 2); 
      //Radiative factor
      double xs_p = grXS_pair->Eval(e_gamma * 1e3);
      double xs_t = grXS_trip->Eval(e_gamma * 1e3);
      if (e_gamma > 100.0) {
	xs_p = grXS_pair->Eval(100.0 * 1e3);
	xs_t = grXS_trip->Eval(100.0 * 1e3);
      }
      double RadiativeFactorConstant = 0.0093;
      double xs_ratio = xs_t / xs_p;
      double RadiativeFactor = 1.0 + RadiativeFactorConstant / xs_ratio;
      
      h_lgam2->Fill(e_gamma, 1.0 / xs);
      
      xs *= (Z * ScreeningFactor * RadiativeFactor); 
      er_xs *= (Z * ScreeningFactor * RadiativeFactor);
      
      h_lgam3->Fill(e_gamma, 1.0 / xs);
      h_Tkin_rec->Fill(e_recoil_4Vec.E());
      h_theta_vs_Tkin_rec->Fill(e_recoil_4Vec.E(), log10(e_recoil_4Vec.Theta() * TMath::RadToDeg()));
      for (int i = 0; i < npart_photon; i ++) {
	h_Tkin_gam->Fill(photon_4Vec[i].E());
	h_theta_vs_Tkin_gam->Fill(photon_4Vec[i].E(), log10(photon_4Vec[i].Theta() * TMath::RadToDeg()));
      }
      
      //HDDM STUFF
      tmpEvt_t tmpEvt;
      tmpEvt.beam = gamma_4Vec;
      tmpEvt.target = Target_4Vec;
      int j = 0;
      for (int i = 2; i < npart; i ++) {
	if (i != nid_recoil) {
	  tmpEvt.q[j] = TLorentzVector(px[i], py[i], pz[i], e[i]);
	  tmpEvt.pdg[j] = pdg[i];
	  j ++;
	}
      }
      tmpEvt.recoil = e_recoil_4Vec;
      tmpEvt.nGen = npart - 2;
      tmpEvt.rxn = m_process;
      tmpEvt.weight = xs;
      hddmWriter->write(tmpEvt, runNum, counter);
    }
  }
  
  h_Tkin_gam->Write();
  h_Tkin_rec->Write();
  h_theta_vs_Tkin_gam->Write();
  h_theta_vs_Tkin_rec->Write();
  h_egam1->Write();
  h_egam2->Write();
  h_lgam1->Write();
  h_lgam2->Write();
  h_lgam3->Write();
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  
  return 0;
}
