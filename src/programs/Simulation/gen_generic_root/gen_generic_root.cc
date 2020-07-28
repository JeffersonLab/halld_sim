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

#include <TGenPhaseSpace.h>

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
  }

  // Get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(genconfigfile);
  Double_t * m_target_mass = ReadFile->GetConfig1Par("target_mass");
  Double_t * m_nb_fs = ReadFile->GetConfig1Par("nb_fs");
  Double_t * m_masses1 = ReadFile->GetConfig6Par("masses1");
  Double_t * m_masses2 = ReadFile->GetConfig6Par("masses2");
  Double_t * m_pdg1 = ReadFile->GetConfig6Par("pdg1");
  Double_t * m_pdg2 = ReadFile->GetConfig6Par("pdg2");
  
  int npart_thrown = (int) m_nb_fs[0];
  double masses[npart_thrown];
  int pdg[npart_thrown];
  double m_threshold = 0;
  if (npart_thrown < 12) {
    for (int i = 0; i < npart_thrown; i ++) {
      if (i < 6) {
	masses[i] = m_masses1[i];
	pdg[i] = m_pdg1[i];
      } if (i >= 6) {
	masses[i] = m_masses2[i - 6];
	pdg[i] = m_pdg2[i - 6];
      }
      m_threshold += masses[i];
    }
  }
  TGenPhaseSpace decayGen;
  TLorentzVector target_vec(0, 0, 0, m_target_mass[0]);
  
  TFile * diagOut = new TFile("gen_generic_root.root", "recreate" );
  TH1F * h_egam = new TH1F("egam", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_wgam = new TH1F("wgam", ";W_{#gammap} [GeV/c^{2}];Count/MeV", (int) ((5.0 - 0.93827200)*1e3), 0.93827200, 5.0);
  
  TH2F * h_theta_vs_Tkin[npart_thrown];
  for (int i = 0; i < npart_thrown; i ++)
    h_theta_vs_Tkin[i] = new TH2F(Form("theta_vs_Tkin_%d", i), Form(";T_{%d}^{kin} [GeV];log_{10}(#theta_{$d}) [^{o}];Count/10MeV", i), 1200, 0.0, 12.0, 1200, -5, 2.25);
  
  for (int i = 0; i < nEvents; i ++) { //Generate photon beam spectrum
    if (i%10000 == 1)
      cout << "event " << i <<endl;
    
    // get beam energy
    double ebeam = 0;
    if (beamconfigfile == "" || cobrem_vs_E == 0) 
      ebeam = ebeam_spectrum.GetRandom();
    else if (beamconfigfile != "")
      ebeam = cobrem_vs_E->GetRandom();
    
    h_egam->Fill(ebeam);
    
    TLorentzVector iphoton_vec = TLorentzVector(0, 0, ebeam, ebeam);
    TLorentzVector is_vec = iphoton_vec + target_vec;
    h_wgam->Fill(is_vec.M());
    
    if (is_vec.M() >= m_threshold) {
      
      //HDDM STUFF
      tmpEvt_t tmpEvt;
      tmpEvt.beam = iphoton_vec;
      tmpEvt.target = target_vec;
      
      TLorentzVector particle_vec[12];
      if (decayGen.SetDecay(is_vec, npart_thrown, masses)) {
	decayGen.Generate();
	for (int j = 0; j < npart_thrown; j ++) {
	  particle_vec[j] = * decayGen.GetDecay(j);
	  tmpEvt.q[j] = particle_vec[j];
	  tmpEvt.pdg[j] = pdg[j];
	  h_theta_vs_Tkin[j]->Fill(particle_vec[j].E() - particle_vec[j].M(), log10(particle_vec[j].Theta() * TMath::RadToDeg()));
	}
      }
      tmpEvt.nGen = npart_thrown;
      tmpEvt.rxn = Form("%d particles phase-space generation", npart_thrown);
      tmpEvt.weight = 1;
      hddmWriter->write(tmpEvt, runNum, i);
    }
  }
  h_egam->Write();
  h_wgam->Write();
  for (int i = 0; i < npart_thrown; i ++)
    h_theta_vs_Tkin[i]->Write();
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  
  return 0;
}
