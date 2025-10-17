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

//#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
//#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"

//#include "AMPTOOLS_AMPS/Compton.h"

//#include "AMPTOOLS_MCGEN/GammaPToXP.h"

//#include "IUAmpTools/AmpToolsInterface.h"
//#include "IUAmpTools/ConfigFileParser.h"

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
#include <TTimeStamp.h>
#include <TGenPhaseSpace.h>

#include "HddmOut.h"
#include "MyReadConfig.h"
#include "HDDM/hddm_s.hpp"

using std::complex;
using namespace std;

#define GAMMA_TYPE 1
#define ELECTRON_TYPE 3

Double_t* m_fixed_ecm;//!

TLorentzVector doCalEnergy(double BeamEnergy, double nucleusMass, double ParticipantMass, double SpectatorMass, TLorentzVector MesonP4, TLorentzVector RecoilP4);
Double_t RelBW(Double_t *x, Double_t *par);

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
  TRandom3* fRandom = new TRandom3();  
  TTimeStamp * time_st = new TTimeStamp();
  double_t timeseed = time_st->GetNanoSec();
  // random number initialization (set to 0 by default)
  fRandom->SetSeed(timeseed);
  
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
  //Double_t * m_target_mass = ReadFile->GetConfig1Par("target_mass");
  Double_t * m_nb_fs = ReadFile->GetConfig1Par("nb_fs");
  Double_t * m_nb_fsc1 = ReadFile->GetConfig1Par("nb_c1");
  Double_t * m_nb_fsc2 = ReadFile->GetConfig1Par("nb_c2");
  Double_t * m_eg_range = ReadFile->GetConfig2Par("eg_range");
  Double_t * m_th_range = ReadFile->GetConfig2Par("th_range");
  Double_t * m_tk_range = ReadFile->GetConfig2Par("tk_range");
  if (beamconfigfile != "") {
    cobrem_vs_E->GetXaxis()->SetRangeUser(m_eg_range[0], m_eg_range[1]);
  }
  TString * s_particles = ReadFile->GetConfig5str("decay");
  TString * s_particles_child1 = ReadFile->GetConfig5str("child1");
  TString * s_particles_child2 = ReadFile->GetConfig5str("child2");
  int npart_thrown_c1 = (int) m_nb_fsc1[0];
  int npart_thrown_c2 = (int) m_nb_fsc2[0];
      
  if (ReadFile->GetConfigName("fixed_ecm") != "") {
    m_fixed_ecm = ReadFile->GetConfig2Par("fixed_ecm");
    cout <<"Low ecm " << m_fixed_ecm[0] << " high ecm " << m_fixed_ecm[1] << endl;
  }

  double masses_c1[npart_thrown_c1];
  double masses_c2[npart_thrown_c2];
  
  int npart_thrown = (int) m_nb_fs[0];
  Particle_t t_particles[npart_thrown + npart_thrown_c1 + npart_thrown_c2];
  Particle_t t_dummy;
  double masses[npart_thrown];
  int pdg[npart_thrown + npart_thrown_c1 + npart_thrown_c2];
  double m_threshold = 0;
  TF1 * fBW = NULL;
  bool b_bw = false;
  cout << "Number of particles thrown " << npart_thrown << endl;
  if (npart_thrown < 5) {
    int j = 0;
    for (int i = 0; i < npart_thrown; i ++) {
      if (!s_particles[i].Contains("980") && !s_particles[i].Contains("1320")) {
	t_particles[j] = ParticleEnum(s_particles[i].Data());
	masses[i] = ParticleMass(t_particles[j]);
	pdg[j] = PDGtype(t_particles[j]);
	m_threshold += masses[i];
	cout <<"index " << i << " s_particles[" << j << "] " << s_particles[i] << " mass " << masses[i] << " pdg " << pdg[j] << endl;
	j ++;
      }
      if (s_particles[i].Contains("a") && s_particles[i].Contains("980")) {
	double mass_sum = 0;
	for (int i = 0; i < npart_thrown_c1; i ++) {
	  t_dummy = ParticleEnum(s_particles_child1[i].Data());
	  mass_sum += ParticleMass(t_dummy);
	}
	fBW = new TF1("fBW", RelBW, mass_sum, 2.0, 3);
	fBW->SetParameters(0.98, 0.075, 1.0);
	b_bw = true;
	cout << "Loading BW TF1 " << s_particles[i] << endl;
      }
      if (s_particles[i].Contains("a") && s_particles[i].Contains("1320")) {
	double mass_sum = 0;
	for (int i = 0; i < npart_thrown_c1; i ++) {
	  t_dummy = ParticleEnum(s_particles_child1[i].Data());
	  mass_sum += ParticleMass(t_dummy);
	}
	fBW = new TF1("fBW", RelBW, mass_sum, 2.5, 3);
	fBW->SetParameters(1.32, 0.105, 1.0);
	b_bw = true;
	cout << "Loading BW TF1 " << s_particles[i] << endl;
      }
    }
    for (int i = 0; i < npart_thrown_c1; i ++) {
      t_particles[j] = ParticleEnum(s_particles_child1[i].Data());
      masses_c1[i] = ParticleMass(t_particles[j]);
      pdg[j] = PDGtype(t_particles[j]);
      cout << "index " << i << " s_particles[" << j << "] " << s_particles_child1[i] << " mass " << masses_c1[i] << " pdg " << pdg[j] << endl;
      j ++;
    }
    for (int i = 0; i < npart_thrown_c2; i ++) {
      t_particles[j] = ParticleEnum(s_particles_child2[i].Data());
      masses_c2[i] = ParticleMass(t_particles[j]);
      pdg[j] = PDGtype(t_particles[j]);
      cout << "index " << i << " s_particles[" << j << "] " << s_particles_child2[i] << " mass " << masses_c2[i] << " pdg " << pdg[j] << endl;
      j ++;
    }
  }
  if (npart_thrown_c1 > 0 && npart_thrown_c2 == 0) {
    for (int i = 0; i < (npart_thrown -1 + npart_thrown_c1); i ++) {
      cout  << "index " << i << " mass[" << "] = " << ParticleMass(t_particles[i]) << " particle " << t_particles[i] << endl;
    }
  } else if (npart_thrown_c1 > 0 && npart_thrown_c2 > 0) {
    for (int i = 0; i < (npart_thrown_c1 + npart_thrown_c2); i ++) {
      cout  << "index " << i << " mass[" << "] = " << ParticleMass(t_particles[i]) << " particle " << t_particles[i] << endl;
    }
  }
  bool b_Free = true;
  bool b_Bound = false;
  TString m_target = "";
  TString m_Fermi_file = "";
  TString m_Participant = "";
  TString m_Spectator = "";
  TH1F * m_h_PFermi = new TH1F("PFermi", "", 12000, 0.0, 12.0);
  TH2F * h_sf = NULL;
  Particle_t t_target;
  Particle_t t_spectator;
  Particle_t t_participant;
  if (ReadFile->GetConfigName("target") != "") {
    m_target = ReadFile->GetConfigName("target");
    cout << "target " << m_target << endl;
  }
  if (ReadFile->GetConfigName("fermi_file") != "" && ReadFile->GetConfigName("participant") != "" && ReadFile->GetConfigName("spectator") != "") {
    b_Free = false;
    b_Bound = true;
    m_Fermi_file = ReadFile->GetConfigName("fermi_file");
    cout << "Fermi_file " << m_Fermi_file << endl;
    m_Participant = ReadFile->GetConfigName("participant");
    m_Spectator = ReadFile->GetConfigName("spectator"); 
    cout <<"Target is made of " << m_target << " with the participant " << m_Participant << " and spectator " << m_Spectator << endl;
    cout <<"Nucleon Fermi motion is located in " << m_Fermi_file << endl;
    if (!m_Fermi_file.Contains("SRC-unweighted")) {
      ifstream in;
      in.open(m_Fermi_file);
      int i = 0;
      while (in.good()) {
	double pf = 0, val = 0;
	in >> pf >> val;
	if (val > 0) {
	  m_h_PFermi->SetBinContent(i + 1, val);
	  i ++;
	}
      }
      in.close();
    } else if (m_Fermi_file.Contains("SRC-unweighted")) {
      TFile * t_sf = new TFile(m_Fermi_file);
      if (m_Participant == "Neutron") h_sf = (TH2F *) t_sf->Get("src_sf_n");
      if (m_Participant == "Proton") h_sf = (TH2F *) t_sf->Get("src_sf_p");
      if (m_Participant == "Deuteron") h_sf = (TH2F *) t_sf->Get("src_sf_d");
    }
    t_target = ParticleEnum(m_target.Data());
    t_spectator = ParticleEnum(m_Spectator.Data());
    t_participant = ParticleEnum(m_Participant.Data());
    cout << "Target mass " << ParticleMass(t_target) << " pdg " << PDGtype(t_target) << endl;
    cout << "Spectator mass " << ParticleMass(t_spectator) << " pdg " << PDGtype(t_spectator) << endl;
    cout << "Participant mass " << ParticleMass(t_participant) << " pdg " << PDGtype(t_participant) << endl;
  } else {
    t_target = ParticleEnum(m_target.Data());
    cout << "Target mass " << ParticleMass(t_target) << " pdg " << PDGtype(t_target) << endl;
  }
  
  TGenPhaseSpace decayGen1;
  TGenPhaseSpace decayGen2;
 
  TFile * diagOut = new TFile("gen_generic_root.root", "recreate" );
  TH1F * h_egam = new TH1F("egam", ";E_{#gamma} [GeV];Count/MeV", 12000, 0.0, 12.0);
  TH1F * h_wgam = new TH1F("wgam", ";W_{#gammap} [GeV/c^{2}];Count/MeV", 11000, 1.0, 12.0);
  TH1F * h_wgam_IS = new TH1F("wgam_IS", ";W_{#gammap} [GeV/c^{2}];Count/MeV", 11000, 1.0, 12.0);
  TH1F * h_wgam_FS = new TH1F("wgam_FS", ";W_{#gammap} [GeV/c^{2}];Count/MeV", 11000, 1.0, 12.0);
  TH1F * h_wgam_FS_check = new TH1F("wgam_FS_check", ";W_{#gammap} [GeV/c^{2}];Count/MeV", 11000, 1.0, 12.0);
  TH1F * h_thrown_FermiP1 = new TH1F("thrown_FermiP1",";p_{F} [GeV/c];",250,0.,1.);
  TH1F * h_thrown_FermiP2 = new TH1F("thrown_FermiP2",";p_{F} [GeV/c];",250,0.,1.);
  TH1F * h_thrown_FermiP3 = new TH1F("thrown_FermiP3",";p_{F} [GeV/c];",250,0.,1.);
  TH1F * h_cons1[2];
  TH1F * h_cons2[2];
  TH1F * h_bw_mass = new TH1F("bw_mass", ";m_{a} [GeV/c^{2}];", 250, 0., 3.);
  TH1F * h_bw_mass_check = new TH1F("bw_mass_check", ";m_{a} [GeV/c^{2}];", 250, 0., 3.);
  for (int i = 0; i < 2; i ++) {
    h_cons1[i] = new TH1F(Form("cons1_%d", i), "#DeltaE/M [GeV];Events",1000,-5., 5.);
    h_cons2[i] = new TH1F(Form("cons2_%d", i), "#DeltaE/M [GeV];Events",1000,-5., 5.);
  }
  TH1F * h_cop1 = new TH1F("cop1", ";|#phi_{#eta}-#phi_{recoil}| [^{o}];Events #", 360, 0., 360.);
  TH1F * h_cop2 = new TH1F("cop2", ";|#phi_{#eta}-#phi_{recoil}| [^{o}];Events #", 360, 0., 360.);
  TH1F * h_cop3 = new TH1F("cop3", ";|#phi_{#eta}-#phi_{recoil}| [^{o}];Events #", 360, 0., 360.);
  TH2F * h_theta_vs_Tkin[npart_thrown + npart_thrown_c1 + npart_thrown_c2];
  for (int i = 0; i < (npart_thrown + npart_thrown_c1 + npart_thrown_c2); i ++)
    //h_theta_vs_Tkin[i] = new TH2F(Form("theta_vs_Tkin_%d", i), Form(";T_{%d}^{kin} [GeV];log_{10}(#theta_{$d}) [^{o}];Count/10MeV", i), 1200, 0.0, 12.0, 1200, -5, 2.25);
    h_theta_vs_Tkin[i] = new TH2F(Form("theta_vs_Tkin_%d", i), Form(";T_{%d}^{kin} [GeV];#theta_{$d} [#circ];Count/10MeV", i), 1200, 0.0, 12.0, 1200, 0., 120.);

  TLorentzVector BeamP4(0, 0, 0, 0);
  TLorentzVector TargetP4(0, 0, 0, 0);
  TLorentzVector ATargetP4(0, 0, 0, ParticleMass(t_target));
  TLorentzVector ISP4(0, 0, 0, 0);
  TLorentzVector FSP4(0, 0, 0, 0);
  TLorentzVector FS_fsP4(0, 0, 0, 0);
  TLorentzVector SpectatorP4(0, 0, 0, 0);
  TLorentzVector ParticipantP4(0, 0, 0, 0);
  TLorentzVector RecoilP4(0, 0, 0, 0);
  double p_Fermi = 0, p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
  //double s_mass = 0, p_mass = 0;
  double SpectatorE = 0, ParticipantE = 0;
  int ctr = 0;
  while (ctr < nEvents) {
  
    //for (int i = 0; i < nEvents; i ++) { //Generate photon beam spectrum
    
    // get beam energy
    double ebeam = 0;
    if (beamconfigfile == "" || cobrem_vs_E == 0) 
      ebeam = ebeam_spectrum.GetRandom();
    else if (beamconfigfile != "")
      ebeam = cobrem_vs_E->GetRandom();
    
    ebeam += fRandom->Uniform(-0.0045, 0.0045);

    if (b_bw) {
      m_threshold = 0;
      for (int j = 0; j < npart_thrown; j ++) {
	if (s_particles[j] == "a0(980)" || s_particles[j] == "a+(980)" || s_particles[j] == "a-(980)" ||
	    s_particles[j] == "a2(1320)" || s_particles[j] == "a+(1320)" || s_particles[j] == "a-(1320)") {
	  masses[j] = fBW->GetRandom();
	  h_bw_mass->Fill(masses[j]);
	  m_threshold += masses[j];
	} else {
	  masses[j] = ParticleMass(t_particles[j]);
	  m_threshold += masses[j];
	}
      }
    }
    if (ebeam < m_eg_range[0] || ebeam >  m_eg_range[1]) continue;
    //cout << "ebeam " << ebeam << endl;
    h_egam->Fill(ebeam);
    //cout << "ebeam " << ebeam << endl;
    BeamP4 = TLorentzVector(0, 0, ebeam, ebeam);
    TargetP4 = TLorentzVector(0, 0, 0, ParticleMass(t_target));
    ISP4 = BeamP4 + TargetP4;
    h_wgam->Fill(ISP4.M());
    FSP4 = TLorentzVector(0, 0, 0, 0);
    FS_fsP4 = TLorentzVector(0, 0, 0, 0);
    if (b_Bound) {
      // IA variables
      double Ebind = 0;
      if (!m_Fermi_file.Contains("SRC-unweighted"))
	p_Fermi = m_h_PFermi->GetRandom();
      else if (m_Fermi_file.Contains("SRC-unweighted"))
	h_sf->GetRandom2(p_Fermi, Ebind);
      h_thrown_FermiP1->Fill(p_Fermi);
      p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
      fRandom->Sphere(p_Fermi_x, p_Fermi_y, p_Fermi_z, p_Fermi);
      if (Ebind == 0.) {
	SpectatorE = sqrt(pow(ParticleMass(t_spectator), 2.0) + pow(p_Fermi ,2.0));
	ParticipantE = ParticleMass(t_target) - SpectatorE;
      } else {
	SpectatorE = sqrt(pow(ParticleMass(t_spectator), 2.0) + pow(p_Fermi ,2.0));
	ParticipantE = ParticleMass(t_participant) - Ebind;
      }
      SpectatorP4 = TLorentzVector(-p_Fermi_x, -p_Fermi_y, -p_Fermi_z, SpectatorE);
      ParticipantP4 = TLorentzVector(p_Fermi_x, p_Fermi_y, p_Fermi_z, ParticipantE);
      ISP4 = BeamP4 + ParticipantP4;
    }
      
    TLorentzVector particle_lab_P4[12];
    TLorentzVector particle_com_P4[12];
    if (ISP4.M() >= m_threshold) {
      //HDDM STUFF
      tmpEvt_t tmpEvt;
      tmpEvt.beam = BeamP4;
      tmpEvt.target = TargetP4;
      tmpEvt.t_target = t_target;
      int k = 0;
      bool b_bad_sim = false;
      if (decayGen1.SetDecay(ISP4, npart_thrown, masses)) {
	decayGen1.Generate();
	for (int j = 0; j < npart_thrown; j ++) {
	  particle_lab_P4[j] = * decayGen1.GetDecay(j);
	  if (std::isnan(particle_lab_P4[j].E())) b_bad_sim = true;
	  FSP4 += particle_lab_P4[j];
	  if (!b_bw) {
	    tmpEvt.q[j] = particle_lab_P4[j];
	    tmpEvt.t_particles[j] = t_particles[j];
	  }
	  if (b_bw && (!s_particles[j].Contains("980") && !s_particles[j].Contains("1320"))) {
	    tmpEvt.q[k] = particle_lab_P4[j];
	    tmpEvt.t_particles[k] = t_particles[k];
	    FS_fsP4 += particle_lab_P4[j];
	    k ++;
	  }
	}
	//if (b_bad_sim) continue;
	if (b_bw && npart_thrown_c1 != 0 && decayGen2.SetDecay(particle_lab_P4[0], npart_thrown_c1, masses_c1)) {
	  decayGen2.Generate();
	  TLorentzVector bwP4(0, 0, 0, 0);
	  b_bad_sim = false;
	  for (int j = 0; j < npart_thrown_c1; j ++) {
	    particle_lab_P4[npart_thrown + j] = * decayGen2.GetDecay(j);
	    if (std::isnan(particle_lab_P4[npart_thrown + j].E())) b_bad_sim = true;
	    tmpEvt.q[k] = particle_lab_P4[npart_thrown + j];
	    tmpEvt.t_particles[k] = t_particles[k];
	    bwP4 += tmpEvt.q[k];
	    FS_fsP4 += particle_lab_P4[npart_thrown + j];
	    k ++;
	  }
	  h_bw_mass_check->Fill(bwP4.M());
	  //if (b_bad_sim) continue;
	}
	//if (b_bw && npart_thrown_c2 != 0 && decayGen2.SetDecay(particle_lab_P4[1], npart_thrown_c2, masses_c2)) {
	//}
	h_cop1->Fill(fabs(particle_lab_P4[0].Phi() - particle_lab_P4[1].Phi()) * TMath::RadToDeg());
	if (b_bw) {
	  h_cop2->Fill(fabs(particle_lab_P4[2].Phi() - particle_lab_P4[1].Phi()) * TMath::RadToDeg());
	  h_cop3->Fill(fabs(particle_lab_P4[3].Phi() - particle_lab_P4[1].Phi()) * TMath::RadToDeg());
	}
	if (b_Bound) {
	  if (!b_bw) {
	    tmpEvt.q[npart_thrown] = SpectatorP4;
	    tmpEvt.t_particles[npart_thrown] = t_spectator;
	  } else if (b_bw) {
	    //cout << "k rank " << k << endl;
	    tmpEvt.q[k] = SpectatorP4;
	    tmpEvt.t_particles[k] = t_spectator;
	  }
	}
	h_cons1[0]->Fill((BeamP4 + ATargetP4).M() - FSP4.M());
	h_cons1[1]->Fill(BeamP4.E() - particle_lab_P4[0].E());
	if (b_bw) {
	  h_cons2[0]->Fill(BeamP4.E() - particle_lab_P4[2].E());
	  h_cons2[1]->Fill(BeamP4.E() - particle_lab_P4[3].E());
	}
      }
      if (b_Free && b_bw) tmpEvt.nGen = k + 1;
      if (b_Free && !b_bw) tmpEvt.nGen = npart_thrown + 1;
      if (b_Bound && b_bw) tmpEvt.nGen = k + 1;
      if (b_Bound && !b_bw) tmpEvt.nGen = npart_thrown + 1;
      
      if (particle_lab_P4[0].Theta() * TMath::RadToDeg() > m_th_range[1]) continue;
      if ((particle_lab_P4[0].E() - particle_lab_P4[0].M()) < m_tk_range[1]) continue;
      
      if (ctr%10000 == 1)
	cout << "event " << ctr << " target mass " << ATargetP4.M() << endl;
      
      for (int j = 0; j < (npart_thrown + npart_thrown_c1 + npart_thrown_c2); j ++) {
	h_theta_vs_Tkin[j]->Fill(particle_lab_P4[j].E() - particle_lab_P4[j].M(), particle_lab_P4[j].Theta() * TMath::RadToDeg());
      }
      h_wgam_IS->Fill((BeamP4 + ATargetP4).M());
      h_wgam_FS->Fill(FSP4.M());
      h_wgam_FS_check->Fill(FS_fsP4.M());
      h_thrown_FermiP2->Fill((FSP4 - (BeamP4 + ATargetP4)).P());
      h_thrown_FermiP3->Fill((FS_fsP4 - (BeamP4 + ATargetP4)).P());
      if (ReadFile->GetConfigName("fixed_ecm") != "") {
	//cout <<"FSP4 mass " << FSP4.M() << " energy " << FSP4.E() << " ebeam " << ebeam << endl;
	if (FSP4.M() < m_fixed_ecm[0] || FSP4.M() > m_fixed_ecm[1]) continue;
      }
      //tmpEvt.rxn = Form("%d particles phase-space generation", npart_thrown);
      tmpEvt.weight = 1;
      if (hddmWriter) hddmWriter->write(tmpEvt, runNum, ctr);
      ctr ++;
    }
  }
  
  h_egam->Write();
  h_wgam->Write();
  h_wgam_IS->Write();
  h_wgam_FS->Write();
  h_wgam_FS_check->Write();
  h_thrown_FermiP1->Write();
  h_thrown_FermiP2->Write();
  h_thrown_FermiP3->Write();
  for (int i = 0; i < 2; i ++) {
    h_cons1[i]->Write();
    h_cons2[i]->Write();
  }
  h_cop1->Write();
  h_cop2->Write();
  h_cop3->Write();
  h_bw_mass->Write();
  h_bw_mass_check->Write();
  for (int i = 0; i < npart_thrown + npart_thrown_c1 + npart_thrown_c2; i ++)
    h_theta_vs_Tkin[i]->Write();
  diagOut->Close();
    
  if (hddmWriter) delete hddmWriter;
  
  return 0;
}
TLorentzVector doCalEnergy(double BeamEnergy, double nucleusMass, double ParticipantMass, double SpectatorMass, TLorentzVector MesonP4, TLorentzVector Recoil_LAB_P4) {
  double E_MesonP4 = MesonP4.E();
  double p_MesonP4_x = MesonP4.X(); 
  double p_MesonP4_y = MesonP4.Y(); 
  double p_MesonP4_z = MesonP4.Z(); 
  double p_MesonP4 = MesonP4.P();
  double phi = Recoil_LAB_P4.Phi();
  double theta = Recoil_LAB_P4.Theta();
  double b = 2.0 * (p_MesonP4_x * cos(phi) * sin(theta) + p_MesonP4_y * sin(phi) * sin(theta) + p_MesonP4_z * cos(theta) - BeamEnergy * cos(theta));
  double c = p_MesonP4 * p_MesonP4 + BeamEnergy * BeamEnergy - 2.0 * BeamEnergy * p_MesonP4_z;
  double d = BeamEnergy + nucleusMass - E_MesonP4;
  double e = pow(SpectatorMass, 2.0) - pow(ParticipantMass, 2.0) - pow(d, 2.0) + c;
  double Delta = 16.0 * pow(d, 2.0) * (pow(e, 2.0) + pow(b * ParticipantMass, 2.0) - pow(d * ParticipantMass * 2.0,2.0));
  TLorentzVector NewParticle(0.0,0.0,0.0,0.0);
  if(Delta>0.) {
    //double sol1 = (2.0 * e * b - sqrt(Delta)) / (2.0 * (4.0 * TMath::Power(d,2.0) - TMath::Power(b,2.0)));
    double sol2 = (2.0 * e * b + sqrt(Delta)) / (2.0 * (4.0 * pow(d, 2.0) - pow(b, 2.0)));
    //cout << "sol1 " << sol1 << endl;
    //cout << "sol2 " << sol2 << endl;
    //if (sol2 <= 0) sol2 = sol1;
    double newpxcal = sol2 * cos(phi) * sin(theta);
    double newpycal = sol2 * sin(phi) * sin(theta);
    double newpzcal = sol2 * cos(theta);
    double energy = sqrt(pow(sol2, 2.0) + pow(ParticipantMass, 2.0));
    NewParticle = TLorentzVector(newpxcal, newpycal, newpzcal, energy); 
  }
  return NewParticle;
}
Double_t RelBW(Double_t *x, Double_t *par) {
  // par[0] = M0 (resonance mass)
  // par[1] = Gamma (width)
  // par[2] = Normalization
  Double_t M   = x[0];
  Double_t M0  = par[0];
  Double_t G   = par[1];
  Double_t norm = par[2];
  
  Double_t numerator   = M0 * G;
  Double_t denominator = (M*M - M0*M0)*(M*M - M0*M0) + (M0*G)*(M0*G);
  
  return norm * numerator / denominator;
}

