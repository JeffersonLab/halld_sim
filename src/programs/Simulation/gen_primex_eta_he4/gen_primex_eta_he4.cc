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
#include <TGenPhaseSpace.h>
#include "TRandom3.h"
#include "TSystem.h"
#include <TTimeStamp.h>
#include "HddmOut.h"
#include "UTILITIES/MyReadConfig.h"
#include "nucleus/nucleus.h"

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

TLorentzVector meson_lab(TLorentzVector ISP4, double m2, double m3, double ThetaLAB, double PhiLAB);
TLorentzVector meson_com_pf(double th, double ph, double E, double p, double sintheta, double costheta);
TLorentzVector doCalEnergy(double BeamEnergy, double nucleusMass, double ParticipantMass, double SpectatorMass, TLorentzVector MesonP4, TLorentzVector RecoilP4);
Double_t RelBW(Double_t *x, Double_t *par);
Double_t RelBW_massdep(Double_t *x, Double_t *par);
double q_pi(double M, double mPi);
//
/*
// Källén function
inline double lambda(double x, double y, double z){
    return x*x + y*y + z*z - 2*(x*y + x*z + y*z);
}
// Build an orthonormal CM basis {x*, y*, z*} where z* is along k_cm.
// If the initial-state plane is ill-defined, pick a safe y*.
inline void BuildCMBasis(const TLorentzVector& k_cm, const TLorentzVector& p_cm,TVector3& xhat, TVector3& yhat, TVector3& zhat)
{
    zhat = k_cm.Vect().Unit();                      // z* along photon in CM
    TVector3 n = k_cm.Vect().Cross(p_cm.Vect());    // normal to initial plane
    if (n.Mag2() < 1e-20) {                         // nearly collinear fallback
        // pick a vector not parallel to zhat
        TVector3 tmp(1,0,0);
        if (std::fabs(zhat.Dot(tmp)) > 0.9) tmp = TVector3(0,1,0);
        yhat = (tmp - tmp.Dot(zhat)*zhat).Unit();
    } else {
        yhat = n.Unit();                             // y* normal to plane
    }
    xhat = yhat.Cross(zhat).Unit();                  // x* completes right-handed basis
}
*/
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

  TF1 * fBW = NULL;
  
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
  TRandom3* fRandom = new TRandom3();  
  TTimeStamp * time_st = new TTimeStamp();
  double_t timeseed = time_st->GetNanoSec();
  // random number initialization (set to 0 by default)
  fRandom->SetSeed(timeseed);
  nucleus * myNucleus = new nucleus(fRandom); 
      
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
  TString m_sc[2] = {"", ""};
  TString m_target = "";
  TString m_Fermi_file = "";
  TString m_Participant = "";
  TString m_Spectator = "";
  Bool_t do_flat_coh = false;
  Bool_t do_flat_qf = false;
  if (ReadFile->GetConfigName("flat_coh") != "") {
    do_flat_coh = true;
    m_flat_coh_angle_min_max_cut = ReadFile->GetConfig2Par("flat_coh");
    cout << "Swicth to flat coherent" << endl;
  }
  //else if (m_target == "") {
  //cout <<"Add something to produce flat coherent" << endl;
  //exit(1);
  //}
  if (ReadFile->GetConfigName("flat_qf") != "") {
    do_flat_qf = true;
    m_flat_coh_angle_min_max_cut = ReadFile->GetConfig2Par("flat_qf");
    cout << "Swicth to flat quasi-free" << endl;
  }
  //else if (m_target == "") {
  //cout <<"Add something to produce flat qfn" << endl;
  //exit(1);
  //}
  if (ReadFile->GetConfigName("rfile") != "") {
    m_rfile = ReadFile->GetConfigName("rfile"); 
    cout << "rfile " << m_rfile << endl;
  } else {
    cout <<"Please add full to root file w/ xs"<<endl; 
    if (!do_flat_coh && !do_flat_qf) exit(1);
  }
  if (ReadFile->GetConfigName("histo") != "") {
    m_histo = ReadFile->GetConfigName("histo"); 
    cout << "histo " << m_histo << endl;
  } else if (m_histo == "") {
    cout <<"Please add histo name" << endl;
    if (!do_flat_coh && !do_flat_qf) exit(1);
  }
  if (ReadFile->GetConfigName("meson") != "") {
    m_meson = ReadFile->GetConfigName("meson");
    cout << "meson " << m_meson << endl;
    if (m_meson == "pi0") m_meson = "Pi0";
    if (m_meson == "eta") m_meson = "Eta";
    if (m_meson == "eta'") m_meson = "EtaPrime";
    if (m_meson == "omega") m_meson = "Omega";
    if (m_meson == "rho0" || m_meson == "Rho0") {
      m_meson = "Rho0";
      cout << " Rho0 special " << endl;
      if (ReadFile->GetConfigName("sc1") != "" && ReadFile->GetConfigName("sc2") != "") {
	cout << " that decays " << endl;
	m_sc[0] = ReadFile->GetConfigName("sc1");
	m_sc[1] = ReadFile->GetConfigName("sc2");
      }
    }
    if (m_meson != "Eta" && m_meson != "EtaPrime" && m_meson != "Pi0" && m_meson != "Omega" && m_meson != "Rho0") {
      cout <<"Wrong meson choice, please choose between Eta, EtaPrime, or Pi0"<<endl;
      exit(1);
    }
  } else if (m_meson == "") {
    cout <<"Please choose between Eta, EtaPrime, Rho0, Omega, or Pi0"<<endl;
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

  
  
  TH1F * m_h_PFermi = new TH1F("PFermi", "", 12000, 0.0, 12.0);
  Particle_t t_target = Gamma;
  Particle_t t_meson = Gamma;
  Particle_t t_sc1 = Gamma;
  Particle_t t_sc2 = Gamma;
  Particle_t t_spectator = Gamma;
  Particle_t t_participant = Gamma;

  if (m_meson == "Rho0") {
    if (m_sc[0] != "" && m_sc[1] != "") {
      cout << m_meson << " decay into " << m_sc[0] << " and " << m_sc[1] << endl;
      t_sc1 = ParticleEnum(m_sc[0].Data());
      t_sc2 = ParticleEnum(m_sc[1].Data());
      cout << "decaying particle 1 " <<  ParticleMass(t_sc1) << " pdg " << PDGtype(t_sc1) << endl;
      cout << "decaying particle 2 " <<  ParticleMass(t_sc2) << " pdg " << PDGtype(t_sc2) << endl;
    }
  }
    
  TH2F * h_sf = NULL;
  if (ReadFile->GetConfigName("fermi_file") != "" && ReadFile->GetConfigName("participant") != "" && ReadFile->GetConfigName("spectator") != "") {
    m_Fermi_file = ReadFile->GetConfigName("fermi_file");
    cout << "Fermi_file " << m_Fermi_file << endl;
    m_Participant = ReadFile->GetConfigName("participant");
    m_Spectator = ReadFile->GetConfigName("spectator"); 
    cout <<"Target is made of " << m_target << " with the participant " << m_Participant << " and spectator " << m_Spectator << endl;
    cout <<"Nucleon Fermi motion is located in " << m_Fermi_file << endl;
    ifstream in;
    if (!m_Fermi_file.Contains("SRC")) {
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
    }
    if (m_Fermi_file.Contains("SRC-unweighted")) {
      cout << m_Fermi_file << endl;
      TFile * t_sf = new TFile(m_Fermi_file);
      if (m_Participant == "Neutron") h_sf = (TH2F *) t_sf->Get("src_sf_n");
      if (m_Participant == "Proton") h_sf = (TH2F *) t_sf->Get("src_sf_p");
      if (m_Participant == "Deuteron") h_sf = (TH2F *) t_sf->Get("src_sf_d");
    }
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

  
  if (m_meson == "Omega") {
    double Gamma = 0.00849;
    double mass = ParticleMass(t_meson);
    fBW = new TF1("fBW", RelBW, -5 * Gamma + mass, +5 * Gamma + mass, 3);
    //fBW = new TF1("fBW", RelBW, 0.4, 1.1, 3);
    fBW->SetParameters(mass, Gamma, 1.0);
  }
  if (m_meson == "Rho0") {
    double Gamma = 0.1474;
    double mass = ParticleMass(t_meson);
    fBW = new TF1("fBW", RelBW, 0.4, 1.1, 3);
    fBW->SetParameters(mass, Gamma, 1.0);
  }

  
  //double M_target = ParticleMass(t_target);
  
  TLorentzVector ATargetP4(0, 0, 0, ParticleMass(t_target));
  TLorentzVector NTargetP4(0, 0, 0, ParticleMass(t_participant));
  
  // Load eta-meson differential cross-section based on Ilya Larin's calculation, see the *.F program in this directory 
  TFile * ifile;
  TH2F * h_dxs = new TH2F();
  if (!do_flat_coh && !do_flat_qf) { 
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
    if (do_flat_coh) cout <<"Flat coherent bin_theta " << bin_theta << " theta_min " << theta_min << " theta_max " << theta_max << endl;
    if (do_flat_qf) cout <<"Flat qf bin_theta " << bin_theta << " theta_min " << theta_min << " theta_max " << theta_max << endl;
  }

  //double M_meson = ParticleMass(t_meson);
  
  TFile* diagOut = new TFile( "gen_primex_eta_he4_diagnostic.root", "recreate" );
  TH2F* h_Tkin_eta_vs_egam = new TH2F("Tkin_eta_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#eta} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_photon_vs_egam = new TH2F("Tkin_photon_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_recoilA_vs_egam = new TH2F("Tkin_recoilA_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 1.0);
  TH2F* h_theta_eta_vs_egam = new TH2F("theta_eta_vs_egam", ";E_{#gamma} [GeV];#theta_{#eta} [GeV];Count [a.u.]", bin_egam, egam_min, egam_max, bin_theta, theta_min, theta_max);
  //TH2F* h_theta_photon_vs_egam = new TH2F("theta_photon_vs_egam", ";E_{#gamma} [GeV];#theta_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  TH2F* h_theta_recoilA_vs_egam = new TH2F("theta_recoilA_vs_egam", ";E_{#gamma} [GeV];#theta_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  TH1F *thrown_FermiP1 = new TH1F("thrown_FermiP1",";p_{F} [GeV/c];",250,0.,1.);;
  TH1F *thrown_FermiP2 = new TH1F("thrown_FermiP2",";p_{F} [GeV/c];",250,0.,1.);;
  TH1F *thrown_FermiP3 = new TH1F("thrown_FermiP3",";p_{F} [GeV/c];",250,0.,1.);;
  
  TH1F * h_cop = new TH1F("cop", ";|#phi_{#eta}-#phi_{recoil}| [^{o}];Events #", 360, 0., 360.);

  TH1F * h_mass_diff = new TH1F("mass_diff", ";#DeltaM [GeV];Events #", 2000, -1., 1.);

  TH1F * h_meson_mass = new TH1F("meson_mass","; Mass [GeV];Events #", 2000, 0., 2.);
  TH1F * h_meson_theta = new TH1F("meson_theta","; cos;Events #", 2000, -1., 1.);
  
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
    TLorentzVector BeamP4(0, 0, ebeam, ebeam);
    
    ATargetP4 = TLorentzVector(0, 0, 0, ParticleMass(t_target));
    NTargetP4 = TLorentzVector(0, 0, 0, ParticleMass(t_participant));
    TLorentzVector ISP4 = BeamP4 + ATargetP4;
    if (do_flat_qf || m_rfile.Contains("free"))
      ISP4 = BeamP4 + NTargetP4;
    
    // Mass in the centre-of-mass frame
    //bool good_evt = true;
    double ThetaLAB = 0;
    if (!do_flat_coh && !do_flat_qf) {
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
      ThetaLAB = fRandom->Uniform(theta_min, theta_max);
      ThetaLAB *= TMath::DegToRad();
    }
    double mass = ParticleMass(t_meson);
    if (m_meson == "Omega" || m_meson == "Rho0") {
      mass = fBW->GetRandom();
    }
    h_meson_mass->Fill(mass);
    TLorentzVector sc1P4_lab(0, 0, 0, 0);
    TLorentzVector sc2P4_lab(0, 0, 0, 0);
    double PhiLAB = fRandom->Uniform(-TMath::Pi(), TMath::Pi());
    TLorentzVector eta_LAB_P4 = meson_lab(ISP4, mass,  ParticleMass(t_target), ThetaLAB, PhiLAB);
    if (do_flat_qf || m_rfile.Contains("free"))
      eta_LAB_P4 = meson_lab(ISP4, mass,  ParticleMass(t_participant), ThetaLAB, PhiLAB);
    TLorentzVector Recoil_LAB_P4 = ISP4 - eta_LAB_P4;
    TLorentzVector eta_COM_P4 = eta_LAB_P4;
    TLorentzVector Recoil_COM_P4 = Recoil_LAB_P4;
    eta_COM_P4.Boost(-ISP4.BoostVector());
    Recoil_COM_P4.Boost(-ISP4.BoostVector());
    double theta_com = eta_COM_P4.Theta();
    double phi_com = eta_COM_P4.Phi();

    if (std::isnan(eta_LAB_P4.E())) {
      cout <<"Initial -nan / what is the eta polar lab. angle " <<  ThetaLAB * TMath::RadToDeg() << endl;
    }
    
    // IA variables
    double p_Fermi = 0, p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
    TLorentzVector SpectatorP4(0, 0, 0, 0);
    TLorentzVector ParticipantP4(0, 0, 0, 0);
    TVector3 FermiP3(0, 0, 0);
    double weight = 1., Ebind = 0;
    double SpectatorE = 0, ParticipantE = 0;
    if (m_Spectator != "") {
      //double s_mass = ParticleMass(t_spectator);
      //double p_mass = ParticleMass(t_participant);
      if (!m_Fermi_file.Contains("SRC")) {
	p_Fermi = m_h_PFermi->GetRandom();
	thrown_FermiP1->Fill(p_Fermi);
	p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
	fRandom->Sphere(p_Fermi_x, p_Fermi_y, p_Fermi_z, p_Fermi);
	SpectatorE = sqrt(pow(ParticleMass(t_spectator), 2.0) + pow(p_Fermi ,2.0));
	ParticipantE = ParticleMass(t_target) - SpectatorE;
	SpectatorP4 = TLorentzVector(-p_Fermi_x, -p_Fermi_y, -p_Fermi_z, SpectatorE);
	ParticipantP4 = TLorentzVector(p_Fermi_x, p_Fermi_y, p_Fermi_z, ParticipantE);
      } else if (m_Fermi_file.Contains("SRC")) {
	if (!m_Fermi_file.Contains("SRC-unweighted")) {
	  if (m_Participant == "Neutron") {
	    do {
	      weight = 1.;
	      if (m_target == "Deuteron") myNucleus->sample_SF_deut_n(weight, p_Fermi, Ebind);
	      if (m_target == "Helium") myNucleus->sample_SF_He_n(weight, p_Fermi, Ebind);
	      if (m_target == "Carbon") myNucleus->sample_SF_C12_n(weight, p_Fermi, Ebind);
	    } while (weight == 0.);
	  }
	  if (m_Participant == "Proton") {
	    do {
	      weight = 1.;
	      if (m_target == "Deuteron") myNucleus->sample_SF_deut_p(weight, p_Fermi, Ebind);
	      if (m_target == "Helium") myNucleus->sample_SF_He_p(weight, p_Fermi, Ebind);
	      if (m_target == "Carbon") myNucleus->sample_SF_C12_p(weight, p_Fermi, Ebind);
	    } while (weight == 0.);
	  }
	  if (m_Participant == "Deuteron") {
	    do {
	      weight = 1.;
	      if (m_target == "Helium") myNucleus->sample_SF_He_d(weight, p_Fermi, Ebind);
	      if (m_target == "Carbon") myNucleus->sample_SF_C12_d(weight, p_Fermi, Ebind);
	    } while (weight == 0.);
	  }
	} else if (m_Fermi_file.Contains("SRC-unweighted")) {
	  h_sf->GetRandom2(p_Fermi, Ebind, fRandom);
	}
	double phi_src = 2. * TMath::Pi() * fRandom->Uniform();
	double cosTheta_src = -1 + 2. * fRandom->Uniform();
	double theta_src = acos(cosTheta_src);
	FermiP3.SetMagThetaPhi(p_Fermi , theta_src, phi_src);
	ParticipantE = ParticleMass(t_participant) - Ebind;
	if (ParticipantE < 0) cout <<"participant negative mass"<<endl;
	//SpectatorE = ParticleMass(t_target) - ParticipantE;
	//SpectatorE = ParticleMass(t_spectator) - Ebind;
	SpectatorE = sqrt(pow(ParticleMass(t_spectator), 2.0) + pow(p_Fermi ,2.0));
	if (SpectatorE < 0) cout <<"spectator negative mass"<<endl;
	ParticipantP4 = TLorentzVector(FermiP3, ParticipantE);
	SpectatorP4 = TLorentzVector(-FermiP3, SpectatorE);
	//SpectatorP4 = TLorentzVector(-FermiP3.X(), -FermiP3.Y(), -FermiP3.Z(), SpectatorE);
	thrown_FermiP1->Fill(p_Fermi, weight);
	//p_Fermi_x = FermiP3.X();
	//p_Fermi_y = FermiP3.Y();
	//p_Fermi_z = FermiP3.Z();
      }
      /*
      double p_cm_x = p_Fermi_x;
      double p_cm_y = p_Fermi_y;
      double p_cm_z = p_Fermi_z + ebeam;
      double p_cm = sqrt(pow(p_cm_x, 2.) + pow(p_cm_y, 2.) + pow(p_cm_z, 2.));
      double sintheta = p_cm_x / p_cm;
      double costheta = p_cm_z / p_cm;
      */    
      if (m_rfile.Contains("free") || do_flat_qf) {
	ISP4 = BeamP4 + ParticipantP4;
	TLorentzVector BeamP4_cm = BeamP4;
	TLorentzVector ParticipantP4_cm = ParticipantP4;
	BeamP4_cm.Boost(-ISP4.BoostVector());
	ParticipantP4_cm.Boost(-ISP4.BoostVector());
	// Rotate to scattering along z-axis
	double rot_phi = BeamP4_cm.Vect().Phi();
	double rot_theta = BeamP4_cm.Vect().Theta();
	double s = ISP4.M();
	if (s < (mass + ParticleMass(t_participant))) continue;
	double E_eta_com = (pow(s, 2.) - pow(ParticleMass(t_participant), 2.) + pow(mass,2.)) / (2. * s);
	double p_eta_com = sqrt(pow(E_eta_com, 2.) - pow(mass, 2.));
	double E_recoil_com = sqrt(pow(p_eta_com, 2.) + pow(ParticleMass(t_participant), 2.));
	//eta_COM_P4 = meson_com_pf(theta_com, phi_com, E_eta_com, p_eta_com, sintheta, costheta);
	//Recoil_COM_P4 = TLorentzVector(- eta_COM_P4.Px(), - eta_COM_P4.Py(), - eta_COM_P4.Pz(), E_recoil_com);
	TVector3 v_cm(0, 0, 0);
	v_cm.SetMagThetaPhi(p_eta_com, theta_com, phi_com);
	v_cm.RotateY(rot_theta);
	v_cm.RotateZ(rot_phi);
	eta_COM_P4 = TLorentzVector(v_cm, E_eta_com);
	Recoil_COM_P4 = TLorentzVector(-v_cm, E_recoil_com);
	
	eta_LAB_P4 = eta_COM_P4;
	Recoil_LAB_P4 = Recoil_COM_P4;
	eta_LAB_P4.Boost(ISP4.BoostVector());
	Recoil_LAB_P4.Boost(ISP4.BoostVector());
	TVector3 p3_recoil = Recoil_LAB_P4.Vect();
	double E_recoil = sqrt(p3_recoil.Mag() * p3_recoil.Mag() + ParticleMass(t_participant) * ParticleMass(t_participant));
	Recoil_LAB_P4 = TLorentzVector(p3_recoil, E_recoil);
	TVector3 p3_eta = eta_LAB_P4.Vect();
	double E_eta = sqrt(p3_eta.Mag() * p3_eta.Mag() + mass * mass);
	eta_LAB_P4 = TLorentzVector(p3_eta, E_eta);
	//TVector3 p3_spe = SpectatorP4.Vect();
	//double E_spe = sqrt(p3_spe.Mag() * p3_spe.Mag() + ParticleMass(t_spectator) * ParticleMass(t_spectator));
	//SpectatorP4 = TLorentzVector(p3_spe, E_spe);
	h_cop->Fill(fabs(eta_LAB_P4.Phi() - Recoil_LAB_P4.Phi()) * TMath::RadToDeg(), weight);

	TLorentzVector meson_LAB_P4 = meson_lab((BeamP4 + NTargetP4), ParticleMass(t_meson),  ParticleMass(t_participant), eta_LAB_P4.Theta(), eta_LAB_P4.Phi());
	h_mass_diff->Fill(eta_LAB_P4.E() - meson_LAB_P4.E(), weight);
	
	if (fabs(Recoil_LAB_P4.M() - ParticleMass(t_participant)) > 1e-13) {
	  cout << "evt nb " << i << " Participant Mass difference " << fabs(Recoil_LAB_P4.M() - ParticleMass(t_participant))
	       << " kinetic energy " << Recoil_LAB_P4.E() - Recoil_LAB_P4.M() 
	       << " mass evtgen " << Recoil_LAB_P4.M() 
	       << " mass partic " << ParticleMass(t_participant) << endl;
	}      
	if ((eta_LAB_P4.M() - mass) > 1e-13) {
	  cout << "evt nb " << i << " Meson Mass difference " << fabs(eta_LAB_P4.M() - ParticleMass(t_meson))
	       << " kinetic energy " << eta_LAB_P4.E() - eta_LAB_P4.M() 
	       << " mass evtgen " << eta_LAB_P4.M() 
	       << " mass partic " << mass << endl;
	} 
	if (fabs(SpectatorP4.M() - ParticleMass(t_spectator)) > 1e-13) {
	  cout << "evt nb " << i << " Spectator Mass difference " << fabs(SpectatorP4.M() - ParticleMass(t_spectator))
	       << " kinetic energy " << SpectatorP4.E() - SpectatorP4.M() 
	       << " mass evtgen " << SpectatorP4.M() 
	       << " mass partic " << ParticleMass(t_spectator) << endl;
	}
	if (std::isnan(eta_LAB_P4.E())) {
	  cout << "After the fold -nan / what is the Fermi motion " << p_Fermi << " what is the polar com angle " << theta_com * TMath::RadToDeg() << endl;
	  //cout << "sintheta " << sintheta << " costheta " << costheta << " s " << s << endl;
	  continue;
	}
	double tkin_spectator = SpectatorP4.E() - SpectatorP4.M();
	if (tkin_spectator <= 0) cout << "Spectator not moving" << endl;
	tkin_spectator = SpectatorP4.E() - ParticleMass(t_spectator);
	if (tkin_spectator <= 0) cout << "Spectator not moving spec" << endl;
      }
    }
        
    if (m_meson == "Rho0") {
      double m_sc1 = ParticleMass(t_sc1); // GeV
      double m_sc2 = ParticleMass(t_sc2); // GeV
      
      // Momentum magnitude in rho rest frame
      double p = sqrt((pow(mass, 2) - pow(m_sc1 + m_sc2, 2)) * (pow(mass, 2) - pow(m_sc1 - m_sc2, 2))) / (2 * mass);
      
      // Sample cosθ with distribution ~ 1 + cos^2θ
      double costh, phi;
      while (true) {
	costh = fRandom->Uniform(-1., 1.);
	double w = 1 + costh * costh;
	if (fRandom->Uniform(0,2) < w) break;
      }
      phi = fRandom->Uniform(0, 2 * TMath::Pi());
      
      double sinth = sqrt(1 - costh * costh);
      double px = p * sinth * cos(phi);
      double py = p * sinth * sin(phi);
      double pz = p * costh;
      
      h_meson_theta->Fill(costh);
      
      // sc1 4-vector
      TLorentzVector sc1P4_com(px, py, pz, sqrt(pow(p, 2) + pow(m_sc1, 2)));
      
      // sc1 4-vector (back-to-back)
      TLorentzVector sc2P4_com(-px, -py, -pz, sqrt(pow(p, 2) + pow(m_sc2, 2)));
      
      // Combine to check conservation
      sc1P4_lab = sc1P4_com;
      sc2P4_lab = sc2P4_com;
      
      sc1P4_lab.Boost(eta_LAB_P4.BoostVector());
      sc2P4_lab.Boost(eta_LAB_P4.BoostVector());
    }
    
    double tkin_meson = eta_LAB_P4.E() - eta_LAB_P4.M();
    double tkin_recoil = Recoil_LAB_P4.E() - Recoil_LAB_P4.M();
    
    if (tkin_meson <= 0) cout << "Meson not moving" <<endl;
    if (tkin_recoil <= 0) cout << "Recoil not moving" << endl;
    tkin_meson = eta_LAB_P4.E() - mass;
    tkin_recoil = Recoil_LAB_P4.E() - ParticleMass(t_participant);
    if (tkin_meson <= 0) cout << "Meson not moving spec" <<endl;
    if (tkin_recoil <= 0) cout << "Recoil not moving spec" << endl;
        
    int ng_max = 0;
        
    // Deduce by energy and mass conservation the recoil nucleus 4Vec
    TLorentzVector He4_LAB_P4 = ISP4 - eta_LAB_P4;
    if (m_target == "Neutron")
      He4_LAB_P4 = TLorentzVector(He4_LAB_P4.Px(), He4_LAB_P4.Py(), He4_LAB_P4.Pz(), sqrt(pow(M_n, 2) + pow(He4_LAB_P4.P(), 2)));
    if (m_target == "Proton")
      He4_LAB_P4 = TLorentzVector(He4_LAB_P4.Px(), He4_LAB_P4.Py(), He4_LAB_P4.Pz(), sqrt(pow(M_p, 2) + pow(He4_LAB_P4.P(), 2)));

    h_Tkin_recoilA_vs_egam->Fill(ebeam, He4_LAB_P4.E() - He4_LAB_P4.M(), weight);
    h_theta_recoilA_vs_egam->Fill(ebeam, He4_LAB_P4.Theta() * TMath::RadToDeg(), weight);
    h_Tkin_eta_vs_egam->Fill(ebeam, eta_LAB_P4.E() - eta_LAB_P4.M(), weight);
    h_theta_eta_vs_egam->Fill(ebeam, eta_LAB_P4.Theta() * TMath::RadToDeg(), weight);
    if (std::isnan(eta_LAB_P4.E())) {
      cout << "Energy of the eta -nan" << endl;
    } else if (hddmWriter /*&& good_evt*/) {
      // ======= HDDM output =========
      tmpEvt_t tmpEvt;
      tmpEvt.str_decay = "";
      tmpEvt.str_target = m_target;
      tmpEvt.beam = BeamP4;
      tmpEvt.target = ATargetP4;
      tmpEvt.str_meson = m_meson;
      tmpEvt.t_targ = t_target;
      tmpEvt.t_meso = t_meson;
      tmpEvt.q1 = eta_LAB_P4;
      //cout <<"I am here 1 "<<endl;
      if (ng_max == 0 && m_Fermi_file == "" && m_sc[0] == "" && m_sc[1] == "") {
	//cout <<"I am here 2 "<<endl;
	tmpEvt.q2 = He4_LAB_P4;
	tmpEvt.nGen = 2;
      } else if (ng_max == 0 && m_Fermi_file != "" && m_meson != "Rho0") {
	//cout <<"I am here 3 "<<endl;
	//cout <<"part x " << ParticipantP4.X() << " y " << ParticipantP4.Y() << " z " << ParticipantP4.Z() << " e " << ParticipantP4.E() << " m " << ParticipantP4.M() << endl;
	//cout <<"spec x " << SpectatorP4.X() << " y " << SpectatorP4.Y() << " z " << SpectatorP4.Z() << " e " << SpectatorP4.E() << " m " << SpectatorP4.M() << endl;
	TLorentzVector NeutronP4 = doCalEnergy(ebeam, ParticleMass(t_target), ParticipantP4.M(), SpectatorP4.M(), eta_LAB_P4, Recoil_LAB_P4);
	//cout <<"cal p " << NeutronP4.P() << " thrown " << Recoil_LAB_P4.P() << endl;
	//cout <<"cal m " << NeutronP4.M() << " thrown " << Recoil_LAB_P4.M() << " target mass " << ParticleMass(t_target) << endl; 
	thrown_FermiP2->Fill((eta_LAB_P4 + Recoil_LAB_P4 - BeamP4 - ATargetP4).P(), weight);
	thrown_FermiP3->Fill((eta_LAB_P4 + NeutronP4 - BeamP4 - ATargetP4).P(), weight);
	tmpEvt.str_spectator = m_Spectator;
	tmpEvt.str_participant = m_Participant;
	tmpEvt.q2 = Recoil_LAB_P4;
	tmpEvt.q3 = SpectatorP4;
	tmpEvt.t_part = t_participant;
	tmpEvt.t_spec = t_spectator;
	tmpEvt.nGen = 3;
      } else if (ng_max == 0 && m_Fermi_file != "" && m_meson == "Rho0") {
	//cout <<"I am here 3 "<<endl;
	//cout <<"part x " << ParticipantP4.X() << " y " << ParticipantP4.Y() << " z " << ParticipantP4.Z() << " e " << ParticipantP4.E() << " m " << ParticipantP4.M() << endl;
	//cout <<"spec x " << SpectatorP4.X() << " y " << SpectatorP4.Y() << " z " << SpectatorP4.Z() << " e " << SpectatorP4.E() << " m " << SpectatorP4.M() << endl;
	TLorentzVector NeutronP4 = doCalEnergy(ebeam, ParticleMass(t_target), ParticipantP4.M(), SpectatorP4.M(), eta_LAB_P4, Recoil_LAB_P4);
	//cout <<"cal p " << NeutronP4.P() << " thrown " << Recoil_LAB_P4.P() << endl;
	//cout <<"cal m " << NeutronP4.M() << " thrown " << Recoil_LAB_P4.M() << " target mass " << ParticleMass(t_target) << endl; 
	thrown_FermiP2->Fill((eta_LAB_P4 + Recoil_LAB_P4 - BeamP4 - ATargetP4).P(), weight);
	thrown_FermiP3->Fill((eta_LAB_P4 + NeutronP4 - BeamP4 - ATargetP4).P(), weight);
	tmpEvt.str_decay = "decaying";
	tmpEvt.str_spectator = m_Spectator;
	tmpEvt.str_participant = m_Participant;
	tmpEvt.q1 = sc1P4_lab;
	tmpEvt.q2 = sc2P4_lab;
	tmpEvt.t_sc1 = t_sc1;
	tmpEvt.t_sc2 = t_sc2;
	tmpEvt.q3 = Recoil_LAB_P4;
	tmpEvt.q4 = SpectatorP4;
	tmpEvt.t_part = t_participant;
	tmpEvt.t_spec = t_spectator;
	tmpEvt.nGen = 4;
      } else if (ng_max == 0 && m_Fermi_file == "" && m_meson == "Rho0") {
	tmpEvt.str_decay = "decaying";
	tmpEvt.q1 = sc1P4_lab;
	tmpEvt.q2 = sc2P4_lab;
	tmpEvt.t_sc1 = t_sc1;
	tmpEvt.t_sc2 = t_sc2;
	tmpEvt.q3 = He4_LAB_P4;
	tmpEvt.nGen = 3;
      }
      tmpEvt.weight = weight;
      hddmWriter->write(tmpEvt,runNum,i);
    }
    if (asciiWriter) {
      // ======= ASCII output =========
      (*asciiWriter)<<runNum<<" "<<i<<" 3"<<endl;
      // photons from the eta
      //(*asciiWriter)<<"0 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      //(*asciiWriter)<<"   "<<0<<" "<<photon_P4[0].Px()<<" "<<photon_P4[0].Py()<<" "<<photon_P4[0].Pz()<<" "<<photon_P4[0].E()<<endl;			
      //(*asciiWriter)<<"1 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      //(*asciiWriter)<<"   "<<0<<" "<<photon_P4[1].Px()<<" "<<photon_P4[1].Py()<<" "<<photon_P4[1].Pz()<<" "<<photon_P4[1].E()<<endl;			
      // Nucleus recoil
      if (m_target == "He4" || m_target == "Helium") (*asciiWriter)<<"2 "<<Helium_TYPE<<" "<<M_He4<<endl;
      if (m_target == "Be9" || m_target == "Beryllium-9") (*asciiWriter)<<"2 "<<Be9_TYPE<<" "<<M_Be9<<endl;
      if (m_target == "Proton") (*asciiWriter)<<"2 "<<Proton_TYPE<<" "<<M_p<<endl;
      if (m_target == "Neutron") (*asciiWriter)<<"2 "<<Neutron_TYPE<<" "<<M_n<<endl;
      (*asciiWriter)<<"   "<<1<<" "<<He4_LAB_P4.Px()<<" "<<He4_LAB_P4.Py()<<" "<<He4_LAB_P4.Pz()<<" "<<He4_LAB_P4.E()<<endl;			
    }
    
  }
  h_Tkin_eta_vs_egam->Write();
  h_Tkin_photon_vs_egam->Write();
  h_Tkin_recoilA_vs_egam->Write();
  h_theta_eta_vs_egam->Write();
  //h_theta_photon_vs_egam->Write();
  h_theta_recoilA_vs_egam->Write();
  thrown_FermiP1->Write();
  thrown_FermiP2->Write();
  thrown_FermiP3->Write();
  //h_dxs->Write();
  //cobrem_vs_Erange->Write();
  h_mass_diff->Write();
  h_meson_mass->Write();
  h_meson_theta->Write();
  h_cop->Write();
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  if (asciiWriter) delete asciiWriter;
  
  return 0;
}

TLorentzVector meson_com_pf(double th, double ph, double E, double p, double sintheta, double costheta) {
  
  TLorentzVector oldP4(p * sin(th) * cos(ph), p * sin(th) * sin(ph), p * cos(th), E);
  
  return TLorentzVector(+ oldP4.Px() * costheta + oldP4.Pz() * sintheta,
			oldP4.Py(),
			- oldP4.Px() * sintheta + oldP4.Pz() * costheta,
			E);
}

TLorentzVector meson_lab(TLorentzVector ISP4, double m2, double m3, double ThetaLAB, double PhiLAB) {
  double p0 = ISP4.P();
  double E0 = ISP4.E();
  //double m2 = M_meson;
  //double m3 = M_target;
  double E_LAB_eta = 0;
  //double PhiLAB = fRandom->Uniform(-TMath::Pi(), TMath::Pi());
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
  double P_LAB_eta = sqrt(pow(E_LAB_eta, 2) - pow(m2, 2));
  double Px_LAB_eta = P_LAB_eta * sin(ThetaLAB) * cos(PhiLAB);
  double Py_LAB_eta = P_LAB_eta * sin(ThetaLAB) * sin(PhiLAB);
  double Pz_LAB_eta = P_LAB_eta * cos(ThetaLAB);
  
  return TLorentzVector(Px_LAB_eta, Py_LAB_eta, Pz_LAB_eta, E_LAB_eta);
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


// M, M0, Gamma0 in MeV; mPi in MeV
double q_pi(double M, double mPi=0.13957018) {
    if (M <= 2.0*mPi) return 0.0;
    double val = 0.5*M * sqrt(1.0 - 4.0*mPi*mPi/(M*M));
    return val;
}

Double_t RelBW_massdep(Double_t *x, Double_t *par) {
  Double_t M   = x[0];         // MeV
  Double_t M0  = par[0];      // MeV
  Double_t G0  = par[1];      // MeV (width at resonance)
  Double_t norm= par[2];

  // pion mass (use charged or neutral as appropriate)
  const double mPi = 0.13957018; // MeV (pi+)
  double q    = q_pi(M, mPi);
  double q0   = q_pi(M0, mPi);
  double GammaM = 0.0;
  if (q0 > 0.0 && q > 0.0) {
      // P-wave (l=1) -> q^3 behavior and  M0/M factor
      GammaM = G0 * pow(q/q0, 3) * (M0 / M);
  } else {
      GammaM = 0.0;
  }

  double numer = M0 * G0; // keep numerator with on-shell width if you want that convention
  double denom = (M*M - M0*M0)*(M*M - M0*M0) + (M0*GammaM)*(M0*GammaM);
  return norm * numer / denom;
}
