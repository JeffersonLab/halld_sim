/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2019       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *  
*
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
#include <TGenPhaseSpace.h>
#include "TRandom3.h"
#include "TSystem.h"

#include "HddmOut.h"
#include "UTILITIES/MyReadConfig.h"

using std::complex;
using namespace std;

#define eta_TYPE 17
#define etaprime_TYPE 35
#define pi0_TYPE 7
#define gamma_TYPE 1
#define Helium_TYPE 47
#define Be9_TYPE 64
#define Proton_TYPE 14
#define Neutron_TYPE 13

Double_t* m_flat_coh_angle_min_max_cut;//!

int main( int argc, char* argv[] ){
  
  string  beamconfigfile("");
  TString genconfigfile("");// = "gen_config.dat";
  string  outname("");
  string  hddmname("");
  
  ifstream in_coherent;
  ifstream in_incoherent;
  
  double beamLowE   = 3.0;
  //double beamHighE  = 12.0;
  double beamHighE  = 22.011;
  const double M_He4 = 3.727379378;
  const double M_Be9 = 8.39479;
  const double M_p = 0.93827208816;
  const double M_n = 0.93956542052;
  //const double M_eta = 0.54730;
  //const double M_etapr = 0.95778;
  //const double M_pi0 = 0.13497685;
  //const double M_gamma = 0.0;
  
  int runNum = 9001;
  int seed = 0;
  
  int nEvents = 10000;

  //int bin_egam = 9400;
  //double egam_min = 2.5995;
  //double egam_max = 12.0005;
  int bin_egam = 1000;
  double egam_min = 12.011;
  double egam_max = 22.011;
  int bin_theta = 650;
  double theta_min = 0.0;
  double theta_max = 6.5;
  
  //double dummy;
  
  //parse command line:
  for (int i = 1; i < argc; i++) {
    
    string arg(argv[i]);
    
    
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
      cout << "\t -o  <name>\t ASCII file output name" << endl;
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
    cout << "No output specificed:  run gen_compton_simple -h for help" << endl;
    exit(1);
  }
  
  if (genconfigfile == "") {
    cout << "No generator configuration file: run gen_primex_eta_he4 -h for help " << endl;
    exit(1);
  }
  
  // random number initialization (set to 0 by default)
  gRandom->SetSeed(seed);
  
  // initialize HDDM output
  HddmOut *hddmWriter = nullptr;
  if (hddmname != "")
    hddmWriter = new HddmOut(hddmname.c_str());
  
  // initialize ASCII output
  ofstream *asciiWriter = nullptr;
  if (outname != "")
    asciiWriter = new ofstream(outname.c_str());
  
  // Get beam properties from configuration file
  TH1D * cobrem_vs_E = 0;
  if (beamconfigfile != "") {
    BeamProperties beamProp( beamconfigfile );
    cobrem_vs_E = (TH1D*)beamProp.GetFlux();
  }
  
  // Get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(genconfigfile);
  TString m_rfile = "";
  TString m_histo = "";
  TString m_meson = "";
  TString m_target = "";
  TString m_Fermi_file = "";
  TString m_Participant = "";
  TString m_Spectator = "";
  Bool_t do_flat_coh = false;
  if (ReadFile->GetConfigName("flat_coh") != "") {
    do_flat_coh = true;
    m_flat_coh_angle_min_max_cut = ReadFile->GetConfig2Par("flat_coh");
    cout << "Swith to flat coherent" << endl;
  } else if (m_target == "") {
    cout <<"Add something to produce flat coherent" << endl;
    //exit(1);
  }
  if (ReadFile->GetConfigName("rfile") != "") {
    m_rfile = ReadFile->GetConfigName("rfile"); 
    cout << "rfile " << m_rfile << endl;
  } else {
    cout <<"Please add full to root file w/ xs"<<endl; 
    if (!do_flat_coh) exit(1);
  }
  if (ReadFile->GetConfigName("histo") != "") {
    m_histo = ReadFile->GetConfigName("histo"); 
    cout << "histo " << m_histo << endl;
  } else if (m_histo == "") {
    cout <<"Please add histo name" << endl;
    if (!do_flat_coh) exit(1);
  }
  if (ReadFile->GetConfigName("meson") != "") {
    m_meson = ReadFile->GetConfigName("meson");
    cout << "meson " << m_meson << endl;
    if (m_meson == "pi0") m_meson = "Pi0";
    if (m_meson == "eta") m_meson = "Eta";
    if (m_meson == "eta'") m_meson = "EtaPrime";
    if (m_meson != "Eta" && m_meson != "EtaPrime" && m_meson != "Pi0") {
      cout <<"Wrong meson choice, please choose between Eta, EtaPrime, or Pi0"<<endl;
      exit(1);
    }
  } else if (m_meson == "") {
    cout <<"Please choose between Eta, EtaPrime, or Pi0"<<endl;
    exit(1);
  }
  if (ReadFile->GetConfigName("target") != "") {
    m_target = ReadFile->GetConfigName("target");
    cout << "target " << m_target << endl;
  } else if (m_target == "") {
    cout <<"Add target" << endl;
    exit(1);
  }
    
  // Assume a beam energy spectrum of 1/E(gamma)
  TF1 ebeam_spectrum("beam_spectrum","1/x",beamLowE,beamHighE);
  cout << "rfile " << m_rfile << endl;
  cout << "histo " << m_histo << endl;
  cout << "meson " << m_meson << endl;
  cout << "target " << m_target << endl;
  cout << "Fermi_file " << m_Fermi_file << endl;

  TH1F * m_h_PFermi = new TH1F("PFermi", "", 1000, 0.0, 1.0);
  Particle_t t_target;
  Particle_t t_meson;
  Particle_t t_spectator = UnknownParticle;
  Particle_t t_participant = UnknownParticle;
  if (ReadFile->GetConfigName("fermi_file") != "" && ReadFile->GetConfigName("participant") != "" && ReadFile->GetConfigName("spectator") != "") {
    m_Fermi_file = ReadFile->GetConfigName("fermi_file");
    cout << "Fermi_file " << m_Fermi_file << endl;
    m_Participant = ReadFile->GetConfigName("participant");
    m_Spectator = ReadFile->GetConfigName("spectator"); 
    cout <<"Target is made of " << m_target << " with the participant " << m_Participant << " and spectator " << m_Spectator << endl;
    cout <<"Nucleon Fermi motion is located in " << m_Fermi_file << endl;
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
    t_target = ParticleEnum(m_target.Data());
    t_meson = ParticleEnum(m_meson.Data());
    t_spectator = ParticleEnum(m_Spectator.Data());
    t_participant = ParticleEnum(m_Participant.Data());
    cout << "Target mass " << ParticleMass(t_target) << " pdg " << PDGtype(t_target) << endl;
    cout << "Meson mass " << ParticleMass(t_meson) << " pdg " << PDGtype(t_meson) << endl;
    cout << "Spectator mass " << ParticleMass(t_spectator) << " pdg " << PDGtype(t_spectator) << endl;
    cout << "Participant mass " << ParticleMass(t_participant) << " pdg " << PDGtype(t_participant) << endl;
  } else {
    t_target = ParticleEnum(m_target.Data());
    t_meson = ParticleEnum(m_meson.Data());
    cout << "Target mass " << ParticleMass(t_target) << " pdg " << PDGtype(t_target) << endl;
    cout << "Meson mass " << ParticleMass(t_meson) << " pdg " << PDGtype(t_meson) << endl;
  }
  
  double M_target = ParticleMass(t_target);
  
  TLorentzVector Target_4Vec(0, 0, 0, ParticleMass(t_target));
  
  // Load eta-meson differential cross-section based on Ilya Larin's calculation, see the *.F program in this directory 
  TFile * ifile;
  TH2F * h_dxs = new TH2F();
  if (!do_flat_coh) { 
    ifile = new TFile(m_rfile);
    h_dxs = (TH2F *) ifile->Get(m_histo);
    bin_egam = h_dxs->GetNbinsX();
    egam_min = h_dxs->GetXaxis()->GetXmin();
    egam_max = h_dxs->GetXaxis()->GetXmax();                                                                                                                                                   
    bin_theta = h_dxs->GetNbinsY();                                                                                                                                                                      
    theta_min = h_dxs->GetYaxis()->GetXmin();
    theta_max = h_dxs->GetYaxis()->GetXmax();   
  } else {
    bin_egam = 700;
    egam_min = 5.0;
    egam_max = 12.0;
    theta_min = m_flat_coh_angle_min_max_cut[0];
    theta_max = m_flat_coh_angle_min_max_cut[1];
    double s_t = 0.001;
    bin_theta = (int) ((theta_max - theta_min) / s_t);
    cout <<"Flat bin_theta " << bin_theta << " theta_min " << theta_min << " theta_max " << theta_max << endl;
  }

  double M_meson = ParticleMass(t_meson);
  
  TFile* diagOut = new TFile( "gen_primex_eta_he4_diagnostic.root", "recreate" );
  TH2F* h_Tkin_eta_vs_egam = new TH2F("Tkin_eta_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#eta} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_photon_vs_egam = new TH2F("Tkin_photon_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_recoilA_vs_egam = new TH2F("Tkin_recoilA_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 1.0);
  TH2F* h_theta_eta_vs_egam = new TH2F("theta_eta_vs_egam", ";E_{#gamma} [GeV];#theta_{#eta} [GeV];Count [a.u.]", bin_egam, egam_min, egam_max, bin_theta, theta_min, theta_max);
  //TH2F* h_theta_photon_vs_egam = new TH2F("theta_photon_vs_egam", ";E_{#gamma} [GeV];#theta_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  TH2F* h_theta_recoilA_vs_egam = new TH2F("theta_recoilA_vs_egam", ";E_{#gamma} [GeV];#theta_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  TH1F *thrown_FermiP = new TH1F("thrown_FermiP",";p_{F} [GeV/c];",250,0.,1.);;

  for (int i = 0; i < nEvents; ++i) {
    if (i%1000 == 1)
      cout << "event " << i <<endl;
    
    // get beam energy
    double ebeam = 0;
    if (beamconfigfile == "" || cobrem_vs_E == 0) 
      ebeam = ebeam_spectrum.GetRandom();
    else if (beamconfigfile != "")
      ebeam = cobrem_vs_E->GetRandom();
    
    // Incident photon-beam 4Vec
    TLorentzVector InGamma_4Vec(0, 0, ebeam, ebeam);
    
    // Initial state 4Vec
    TLorentzVector IS_4Vec = InGamma_4Vec + Target_4Vec;
    
    // Mass in the centre-of-mass frame
    // double sqrt_s = IS_4Vec.M();
    // double s = pow(sqrt_s, 2);
    
    double ThetaLAB = 0;
    if (!do_flat_coh) {
      // Histo. creation that will store the calculated diff. xs. vs. LAB polar angle
      int ebeam_bin = h_dxs->GetXaxis()->FindBin(ebeam);
      TH1F * h_ThetaLAB = (TH1F *) h_dxs->ProjectionY("h_ThetaLAB", ebeam_bin, ebeam_bin);
      if (h_ThetaLAB->GetEntries() == 0) continue;
      // Generate eta-meson theta in LAB
      ThetaLAB = h_ThetaLAB->GetRandom();
      ThetaLAB *= TMath::DegToRad();
      // deletion
      delete h_ThetaLAB;
    } else {
      // Generate eta-meson theta in LAB
      ThetaLAB = gRandom->Uniform(theta_min, theta_max);
      ThetaLAB *= TMath::DegToRad();
    }
    
    // Generate eta-meson phi in LAB
    double PhiLAB = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
    
    // Calculate eta-meson energy in COM
    //double E_COM_eta = (s - pow(M_He4, 2) + pow(M_eta, 2)) / (2.0 * sqrt_s);
    
    // Calculate eta-meson momentum in COM
    //double P_COM_eta = sqrt(pow(E_COM_eta, 2) - pow(M_eta, 2));
    
    // Calculate eta-meson energy in LAB
    double E_LAB_eta = 0;
    double p0 = IS_4Vec.P();
    double E0 = IS_4Vec.E();
    double m2 = M_meson;
    double m3 = M_target;
    double a = pow(2.0 * p0 * cos(ThetaLAB), 2) - pow(2.0 * E0, 2);
    double b = 4.0 * pow(E0, 3) - pow(2.0 * p0, 2) * E0 + pow(2.0 * m2, 2) * E0 - pow(2.0 * m3, 2) * E0;
    double c = 2.0 * pow(p0 * E0, 2) - 2.0 * pow(m2 * E0, 2) + 2.0 * pow(m2 * p0, 2) + 2.0 * pow(m3 * E0, 2) - 2.0 * pow(m3 * p0, 2) + 
      2.0 * pow(m3 * m2, 2) - pow(E0, 4) - pow(p0, 4) - pow(m2, 4) - pow(m3, 4) - pow(2.0 * m2 * p0 * cos(ThetaLAB), 2);
    double E_LAB_eta0 = (-b + sqrt(pow(b, 2) - 4.0 * a * c)) / (2.0 * a); 
    double E_LAB_eta1 = (-b - sqrt(pow(b, 2) - 4.0 * a * c)) / (2.0 * a); 
    if (E_LAB_eta0 < 0) E_LAB_eta = E_LAB_eta1;
    if (E_LAB_eta1 < 0) E_LAB_eta = E_LAB_eta0;
    if (E_LAB_eta0 > E_LAB_eta1) E_LAB_eta = E_LAB_eta0;
    if (E_LAB_eta1 > E_LAB_eta0) E_LAB_eta = E_LAB_eta1;
    
    // Calculate eta-meson momentun im LAB
    double P_LAB_eta = sqrt(pow(E_LAB_eta, 2) - pow(M_meson, 2));
        
    // Calculate the momentum for each direction
    double Px_LAB_eta = P_LAB_eta * sin(ThetaLAB) * cos(PhiLAB);
    double Py_LAB_eta = P_LAB_eta * sin(ThetaLAB) * sin(PhiLAB);
    double Pz_LAB_eta = P_LAB_eta * cos(ThetaLAB);
    
    // Store the results in TLorentzVector for the eta-meson
    TLorentzVector eta_LAB_4Vec(Px_LAB_eta, Py_LAB_eta, Pz_LAB_eta, E_LAB_eta);
    //TLorentzVector eta_COM_4Vec = eta_LAB_4Vec;
    //eta_COM_4Vec.Boost(-IS_4Vec.BoostVector());
    
    // IA variables
    double p_Fermi = 0, p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
    TLorentzVector SpectatorP4(0, 0, 0, 0);
    TLorentzVector ParticipantP4(0, 0, 0, 0);
    if (m_Spectator != "") {
      double s_mass = ParticleMass(t_spectator);
      double p_mass = ParticleMass(t_participant);
      p_Fermi = m_h_PFermi->GetRandom();
      thrown_FermiP->Fill(p_Fermi);
      p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
      gRandom->Sphere(p_Fermi_x, p_Fermi_y, p_Fermi_z, p_Fermi);
      double SpectatorE = sqrt(pow(p_Fermi, 2) + pow(s_mass, 2));
      SpectatorP4 = TLorentzVector(p_Fermi_x, p_Fermi_y, p_Fermi_z, SpectatorE);
      ParticipantP4 = IS_4Vec - eta_LAB_4Vec - SpectatorP4;
      double w_mass = ParticipantP4.M();
      ParticipantP4 = p_mass / w_mass * ParticipantP4;
    }
        
    int ng_max = 0;
        
    // Deduce by energy and mass conservation the recoil nucleus 4Vec
    TLorentzVector He4_LAB_4Vec = IS_4Vec - eta_LAB_4Vec;
    if (m_target == "Neutron")
      He4_LAB_4Vec = TLorentzVector(He4_LAB_4Vec.Px(), He4_LAB_4Vec.Py(), He4_LAB_4Vec.Pz(), sqrt(pow(M_n, 2) + pow(He4_LAB_4Vec.P(), 2)));
    if (m_target == "Proton")
      He4_LAB_4Vec = TLorentzVector(He4_LAB_4Vec.Px(), He4_LAB_4Vec.Py(), He4_LAB_4Vec.Pz(), sqrt(pow(M_p, 2) + pow(He4_LAB_4Vec.P(), 2)));

    h_Tkin_recoilA_vs_egam->Fill(ebeam, He4_LAB_4Vec.E() - He4_LAB_4Vec.M());
    h_theta_recoilA_vs_egam->Fill(ebeam, He4_LAB_4Vec.Theta() * TMath::RadToDeg());
    h_Tkin_eta_vs_egam->Fill(ebeam, eta_LAB_4Vec.E() - eta_LAB_4Vec.M());
    h_theta_eta_vs_egam->Fill(ebeam, eta_LAB_4Vec.Theta() * TMath::RadToDeg());
            
    if (hddmWriter) {
      // ======= HDDM output =========
      tmpEvt_t tmpEvt;
      tmpEvt.str_target = m_target;
      tmpEvt.beam = InGamma_4Vec;
      tmpEvt.target = Target_4Vec;
      tmpEvt.str_meson = m_meson;
      tmpEvt.t_targ = t_target;
      tmpEvt.t_meso = t_meson;
      tmpEvt.q1 = eta_LAB_4Vec;
      if (ng_max == 0 && m_Fermi_file == "") {
	tmpEvt.q2 = He4_LAB_4Vec;
	tmpEvt.nGen = 2;
      } else if (ng_max == 0 && m_Fermi_file != "") {
	tmpEvt.str_spectator = m_Spectator;
	tmpEvt.str_participant = m_Participant;
	tmpEvt.q2 = ParticipantP4;
	tmpEvt.q3 = SpectatorP4;
	tmpEvt.t_part = t_participant;
	tmpEvt.t_spec = t_spectator;
	tmpEvt.nGen = 3;
      }
      tmpEvt.weight = 1.;
      hddmWriter->write(tmpEvt,runNum,i);
    }
    if (asciiWriter) {
      // ======= ASCII output =========
      (*asciiWriter)<<runNum<<" "<<i<<" 3"<<endl;
      // photons from the eta
      //(*asciiWriter)<<"0 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      //(*asciiWriter)<<"   "<<0<<" "<<photon_4Vec[0].Px()<<" "<<photon_4Vec[0].Py()<<" "<<photon_4Vec[0].Pz()<<" "<<photon_4Vec[0].E()<<endl;			
      //(*asciiWriter)<<"1 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      //(*asciiWriter)<<"   "<<0<<" "<<photon_4Vec[1].Px()<<" "<<photon_4Vec[1].Py()<<" "<<photon_4Vec[1].Pz()<<" "<<photon_4Vec[1].E()<<endl;			
      // Nucleus recoil
      if (m_target == "He4" || m_target == "Helium") (*asciiWriter)<<"2 "<<Helium_TYPE<<" "<<M_He4<<endl;
      if (m_target == "Be9" || m_target == "Beryllium-9") (*asciiWriter)<<"2 "<<Be9_TYPE<<" "<<M_Be9<<endl;
      if (m_target == "Proton") (*asciiWriter)<<"2 "<<Proton_TYPE<<" "<<M_p<<endl;
      if (m_target == "Neutron") (*asciiWriter)<<"2 "<<Neutron_TYPE<<" "<<M_n<<endl;
      (*asciiWriter)<<"   "<<1<<" "<<He4_LAB_4Vec.Px()<<" "<<He4_LAB_4Vec.Py()<<" "<<He4_LAB_4Vec.Pz()<<" "<<He4_LAB_4Vec.E()<<endl;			
    }
    
  }
  h_Tkin_eta_vs_egam->Write();
  h_Tkin_photon_vs_egam->Write();
  h_Tkin_recoilA_vs_egam->Write();
  h_theta_eta_vs_egam->Write();
  //h_theta_photon_vs_egam->Write();
  h_theta_recoilA_vs_egam->Write();
  thrown_FermiP->Write();
  //h_dxs->Write();
  //cobrem_vs_Erange->Write();
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  if (asciiWriter) delete asciiWriter;
  
  return 0;
}


