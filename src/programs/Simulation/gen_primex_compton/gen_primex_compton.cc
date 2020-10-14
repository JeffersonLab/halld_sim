/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       HallD                                           * 
*                                                                         *                                                                                                                               
* Author: The HallD                                                       *                                                                                                                                
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
#include "UTILITIES/MyReadConfig.h"
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
  TLorentzVector Target_4Vec(0, 0, 0, M_electron);
  
  int runNum = 9001;
  int seed = 0;
  
  int nEvents = 10000;
  
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
  
  if (outname.size() == 0 && hddmname == "") {
    cout << "No output specificed:  run gen_whizard -h for help" << endl;
    exit(1);
  }
  
  if (genconfigfile == "") {
    cout << "No generator configuration file: run gen_whizard -h for help " << endl;
    exit(1);
  }
  // random number initialization (set to 0 by default)
  gRandom->SetSeed(seed);
  
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
    cout << cobrem_vs_E->GetEntries() << endl;
  }
  
  // Get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(genconfigfile);
  TString m_run = ReadFile->GetConfigName("run"); 
  TString m_workflow = ReadFile->GetConfigName("workflow"); 
  TString m_file = ReadFile->GetConfigName("file");
  TString m_dir = ReadFile->GetConfigName("dir");
  TString m_out_dir = ReadFile->GetConfigName("out_dir"); 
  TString m_shell = ReadFile->GetConfigName("shell");
  Double_t * m_Ee = ReadFile->GetConfig1Par("Ee");  
  TString m_tagger_file = ReadFile->GetConfigName("tagger_file");

  TFile * diagOut = new TFile( TString::Format("gen_primex_compton_runnb_%d.root", runNum), "recreate" );
  TH1F * h_egam1 = new TH1F("egam1", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_egam2 = new TH1F("egam2", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_Tkin_gam = new TH1F("Tkin_gam", ";T_{#gamma}^{kin} [GeV];Count/10MeV", 1200, 0.0, 12.0);
  TH1F * h_Tkin_rec = new TH1F("Tkin_rec", ";T_{e^{-}-recoil}^{kin} [GeV];Count/10MeV", 1200, 0.0, 12.0);
  TH2F * h_theta_vs_Tkin_gam = new TH2F("theta_vs_Tkin_gam", ";T_{#gamma}^{kin} [GeV];log_{10}(#theta) [^{o}];Count/10MeV", 1200, 0.0, 12.0, 1200, -5, 2.25);
  TH2F * h_theta_vs_Tkin_rec = new TH2F("theta_vs_Tkin_rec", ";T_{e^{-}-recoil}^{kin} [GeV];log_{10}(#theta) [^{o}];Count/10MeV", 1200, 0.0, 12.0, 1200, -5, 2.25);
 
  if (m_run == "true") { //If "true" run WO
    if (m_tagger_file = "") {
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
    } else if (m_tagger_file.Contains("primex_tagm.txt") || m_tagger_file.Contains("primex_tagh.txt")) {
      
      int tagger_channel_nb = 102;
      if (m_tagger_file.Contains("primex_tagh.txt"))
	tagger_channel_nb = 274;

      if (m_workflow != "" && m_shell == "bash")
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.sh %s %f %d %d %d %s %s %d", 
			       //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum, m_workflow.Data()));
			       m_out_dir.Data(), m_Ee[0], nEvents, runNum, runNum, m_workflow.Data(), m_tagger_file.Data(), tagger_channel_nb));
      else if (m_workflow != "" && m_shell == "tcsh")
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.csh %s %f %d %d %d %s %s %d", 
			       //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum, m_workflow.Data()));
			       m_out_dir.Data(), m_Ee[0], nEvents, runNum, runNum, m_workflow.Data(), m_tagger_file.Data(), tagger_channel_nb));
      else if (m_workflow == "" && m_shell == "bash")
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d %s %d", 
			       //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum));
			       m_out_dir.Data(), m_Ee[0], nEvents, runNum, runNum, m_tagger_file.Data(), tagger_channel_nb));
      else if (m_workflow == "" && m_shell == "tcsh")
	system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d %s %d", 
			       //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum));
			       m_out_dir.Data(), m_Ee[0], nEvents, runNum, runNum, m_tagger_file.Data(), tagger_channel_nb));
    }
    
    for (int i = 0; i < h_egam1->GetNbinsX(); i ++) { //Generate LHE file
      double egam = h_egam1->GetBinCenter(i + 1);
      int nbofevt =  h_egam1->GetBinContent(i + 1);
      if (nbofevt > 0) {
	if (m_workflow != "" && m_shell == "bash")
	  system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.sh %s %f %d %d %d %s test 1", 
				 //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum, m_workflow.Data()));
				 m_out_dir.Data(), egam, nbofevt, runNum, runNum, m_workflow.Data()));
	else if (m_workflow != "" && m_shell == "tcsh")
	  system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_slurm.csh %s %f %d %d %d %s test 1", 
				 //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum, m_workflow.Data()));
				 m_out_dir.Data(), egam, nbofevt, runNum, runNum, m_workflow.Data()));
	else if (m_workflow == "" && m_shell == "bash")
	  system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d test 1", 
				 //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum));
				 m_out_dir.Data(), egam, nbofevt, runNum, runNum));
	else if (m_workflow == "" && m_shell == "tcsh")
	  system(TString::Format("source $HALLD_SIM_HOME/src/programs/Simulation/gen_primex_compton/run/compton_prompt.sh %s %f %d %d %d test 1", 
				 //m_out_dir.Data(), (int) egam, nbofevt, runNum, runNum));
				 m_out_dir.Data(), egam, nbofevt, runNum, runNum));
      }
    }
  }
  
  if (m_run == "false" || m_dir != "" || m_file != "") { //Read and loop over a single or all root file
    
    TSystemDirectory dir(m_dir.Data(), m_dir.Data());
    TList * files = dir.GetListOfFiles(); 
    
    float xs = 0;
    float Ebeam = 0;
    float Egamma1 = 0, Egamma2 = 0, Ee = 0;
    float ctheta1 = 0, ctheta2 = 0, costhe = 0;
    float phi1 = 0, phi2 = 0, phie = 0;
    float cswitch = 0, hardcorr = 0, vscsec = 0;

    TChain * m_tree = new TChain("h21");
    
    m_tree->SetBranchAddress("Ebeam", &Ebeam);
    m_tree->SetBranchAddress("Egamma1", &Egamma1);
    m_tree->SetBranchAddress("Egamma2", &Egamma2);
    m_tree->SetBranchAddress("Ee", &Ee);
    m_tree->SetBranchAddress("ctheta1", &ctheta1);
    m_tree->SetBranchAddress("ctheta2", &ctheta2);
    m_tree->SetBranchAddress("costhe", &costhe);
    m_tree->SetBranchAddress("phi1", &phi1);
    m_tree->SetBranchAddress("phi2", &phi2);
    m_tree->SetBranchAddress("phie", &phie);
    m_tree->SetBranchAddress("cswitch", &cswitch);
    m_tree->SetBranchAddress("hardcorr", &hardcorr);
    m_tree->SetBranchAddress("vscsec", &vscsec);

    if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next()) ) {
	fname = file->GetName();
	if (!file->IsDirectory() && fname.Contains(".root")) {
	  TString RootFileName = m_dir + fname;
	  m_tree->Add(RootFileName);
	}
      }
    }
    
    if (m_file != "") 
      m_tree->Add(m_file);
    
    Long64_t NbOfEvent = m_tree->GetEntries();
    TString m_process = "";
    for (int counter = 0; counter < NbOfEvent; counter ++) {
      
      m_tree->GetEntry(counter);
      
      TLorentzVector gamma_4Vec(0, 0, Ebeam, Ebeam);
      TLorentzVector photon_4Vec[10];
      int npart_photon = 1;
      int npart_electron = 1;
      xs = vscsec;
      float ee = Ee;
      float pe = sqrt(pow(ee, 2) - pow(M_electron, 2));
      float thetae = TMath::ACos(costhe);
      float pex = pe * sin(thetae) * cos(phie);
      float pey = pe * sin(thetae) * sin(phie);
      float pez = pe * cos(thetae);
      TLorentzVector e_recoil_4Vec(pex, pey, pez, ee);
      //cout <<"electron px " << pex << " py " << pey << " pz " << pez << " e " << ee << endl; 
      float e1 = Egamma1;
      float e2 = Egamma2;
      float theta1 = TMath::ACos(ctheta1);
      float theta2 = TMath::ACos(ctheta2);
      float p1x = e1 * sin(theta1) * cos(phi1);
      float p1y = e1 * sin(theta1) * sin(phi1);
      float p1z = e1 * cos(theta1);
      float p2x = e2 * sin(theta2) * cos(phi2);
      float p2y = e2 * sin(theta2) * sin(phi2);
      float p2z = e2 * cos(theta2);
      //cout <<"photon1 px " << p1x << " py " << p1y << " pz " << p1z << " e " << e1 << endl; 
      //cout <<"photon1 px " << p2x << " py " << p2y << " pz " << p2z << " e " << e2 << endl; 
      photon_4Vec[0] = TLorentzVector(p1x, p1y, p1z, e1);
      if (e2 > 0) {
	photon_4Vec[1] = TLorentzVector(p2x, p2y, p2z, e2);
	npart_photon = 2;
      }
      TLorentzVector moTransfer = e_recoil_4Vec - Target_4Vec;
      float e_gamma = gamma_4Vec.E();
      h_egam2->Fill(e_gamma);
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
      int npart_thrown = npart_photon + npart_electron;
      for (int i = 0; i < npart_photon; i ++) {
	tmpEvt.q[i] = photon_4Vec[i];
	tmpEvt.pdg[i] = 22;
      }
      if (npart_thrown == 2)
	m_process = "ae_to_ae";
      else if (npart_thrown == 3)
	m_process = "ae_to_aae";
      tmpEvt.recoil = e_recoil_4Vec;
      tmpEvt.nGen = npart_thrown;
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
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  
  return 0;
}
