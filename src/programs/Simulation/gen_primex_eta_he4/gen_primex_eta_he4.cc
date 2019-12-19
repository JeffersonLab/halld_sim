/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2019       GlueX and PrimEX-D Collaborations               * 
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
#include <TGenPhaseSpace.h>
#include "TRandom3.h"
#include "TSystem.h"

#include "HddmOut.h"
#include "MyReadConfig.h"

using std::complex;
using namespace std;

#define eta_TYPE 17
#define gamma_TYPE 22
#define Helium_TYPE 1000020040

int main( int argc, char* argv[] ){
  
  string  beamconfigfile("");
  TString genconfigfile("");// = "gen_config.dat";
  string  outname("");
  string  hddmname("");
  
  ifstream in_coherent;
  ifstream in_incoherent;
  
  double beamLowE   = 3.0;
  double beamHighE  = 12.0;
  const double M_He4 = 3.727379378;
  const double M_eta = 0.54730;
  const double M_gamma = 0.0;
  TLorentzVector Target_4Vec(0, 0, 0, M_He4);
  
  int runNum = 9001;
  int seed = 0;
  
  int nEvents = 10000;
  
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
  
  // initialize HDDM output
  HddmOut *hddmWriter = nullptr;
  if (hddmname != "")
    hddmWriter = new HddmOut(hddmname.c_str());
  
  // initialize ASCII output
  ofstream *asciiWriter = nullptr;
  if (outname != "")
    asciiWriter = new ofstream(outname.c_str());
  
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
  TString m_rfile = ReadFile->GetConfigName("rfile"); 
  TString m_histo = ReadFile->GetConfigName("histo"); 
  
  // Load eta-meson differential cross-section based on Ilya Larin's calculation, see the *.F program in this directory 
  TFile * ifile = new TFile(m_rfile);
  TH2F * h_dxs = (TH2F *) ifile->Get(m_histo);
  
  // Create decayGen
  TGenPhaseSpace decayGen;
  
  TFile* diagOut = new TFile( "gen_primex_eta_he4_diagnostic.root", "recreate" );
  TH2F* h_Tkin_eta_vs_egam = new TH2F("Tkin_eta_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#eta} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_photon_vs_egam = new TH2F("Tkin_photon_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 12.0);
  TH2F* h_Tkin_recoilA_vs_egam = new TH2F("Tkin_recoilA_vs_egam", ";E_{#gamma} [GeV];T^{kin}_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0.0, 1.0);
  TH2F* h_theta_eta_vs_egam = new TH2F("theta_eta_vs_egam", ";E_{#gamma} [GeV];#theta_{#eta} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 10.);
  TH2F* h_theta_photon_vs_egam = new TH2F("theta_photon_vs_egam", ";E_{#gamma} [GeV];#theta_{#gamma} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  TH2F* h_theta_recoilA_vs_egam = new TH2F("theta_recoilA_vs_egam", ";E_{#gamma} [GeV];#theta_{A} [GeV];Count [a.u.]", 1000, 0.0, 12.0, 1000, 0., 180.);
  
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
    double sqrt_s = IS_4Vec.M();
    double s = pow(sqrt_s, 2);
    
    // Histo. creation that will store the calculated diff. xs. vs. LAB polar angle
    int ebeam_bin = h_dxs->GetXaxis()->FindBin(ebeam);
    TH1F * h_ThetaLAB = (TH1F *) h_dxs->ProjectionY("h_ThetaLAB", ebeam_bin, ebeam_bin);
    if (h_ThetaLAB->GetEntries() == 0) continue;
    
    /*
    // XS initialization
    double xs_tot = 0, xs_Primakoff = 0, xs_interference = 0, xs_strong_coherent = 0, xs_incoherent = 0;
    
    // read or calculate differential cross-section
    // open coh. diff. xs
    in_coherent.open(TString::Format("ds_eta_he4_%0.3f.dat", ebeam)); 
    
    // open incoh. diff. xs
    in_incoherent.open(TString::Format("inc_eta_he4_%0.3f.dat", ebeam)); 
    
    // if already calculated, read
    int j = 0; 
    if (in_coherent.good() && in_incoherent.good()) { 
    
    // Read Primakoff, interference, and strong coherent diff. xs terms
    j = 0;
    while (in_coherent.good()) {
    xs_tot = 0; xs_Primakoff = 0; xs_interference = 0; xs_strong_coherent = 0; xs_incoherent = 0;
    in_coherent >> dummy >> dummy >> xs_Primakoff >> xs_interference >> xs_strong_coherent >> xs_incoherent;
    xs_tot = xs_Primakoff + xs_interference + xs_strong_coherent + xs_incoherent;
    h_ThetaLAB->Fill(h_ThetaLAB->GetBinCenter(j + 1), xs_tot);
    j ++;
    }
    in_coherent.close();
    
    // Read incoh. term
    j = 0;
    while (in_incoherent.good()) {
    xs_incoherent = 0;
    in_incoherent >> dummy >> xs_incoherent >> dummy >> dummy;
    h_ThetaLAB->Fill(h_ThetaLAB->GetBinCenter(j + 1), xs_incoherent);
    j ++;
    }
    in_incoherent.close();
    
    } else { // If not calculated, calculate and then read
    
    // Close files if not already closed
    in_coherent.close();
    in_incoherent.close();
    
    // Calculate and produce diff. xs dat files for 100 bin from 0 to 10 degree
    // Calculate Coulomb term of the coherent diff. xs 
    system(TString::Format("ff_coulomb %0.03f ff_coulom_%0.3f.dat", ebeam, ebeam));
    
    // Calculate strong term of the coherent diff. xs 
    system(TString::Format("ff_strong %0.03f ff_strong_%0.3f.dat", ebeam, ebeam));
    
    // Calculate Primakoff, interference, and strong coherent diff. xs 
    system(TString::Format("ds_eta_he4 %0.03f ff_coulom_%0.3f.dat ff_strong_%0.3f.dat ds_eta_he4_%0.3f.dat", ebeam, ebeam, ebeam, ebeam));
    
    // Calculate incoherent
    system(TString::Format("inc_eta_he4 %0.03f inc_eta_he4_%0.3f.dat", ebeam, ebeam));
    
    // Read Primakoff, interference, and strong coherent diff. xs terms
    in_coherent.open(TString::Format("ds_eta_he4_%0.3f.dat", ebeam));
    j = 0;
    while (in_coherent.good()) {
    xs_tot = 0; xs_Primakoff = 0; xs_interference = 0; xs_strong_coherent = 0; xs_incoherent = 0;
    in_coherent >> dummy >> dummy >> xs_Primakoff >> xs_interference >> xs_strong_coherent >> xs_incoherent;
    xs_tot = xs_Primakoff + xs_interference + xs_strong_coherent + xs_incoherent;
    h_ThetaLAB->Fill(h_ThetaLAB->GetBinCenter(j + 1), xs_tot);
    j ++;
    }
    in_coherent.close();
    
    // Read incoh. term
    j = 0;
    while (in_incoherent.good()) {
    xs_incoherent = 0;
    in_incoherent >> dummy >> xs_incoherent >> dummy >> dummy;
    h_ThetaLAB->Fill(h_ThetaLAB->GetBinCenter(j + 1), xs_incoherent);
    j ++;
    }
    in_incoherent.close();
    }
    */
    // Generate eta-meson theta in LAB
    double ThetaLAB = h_ThetaLAB->GetRandom();
    ThetaLAB *= TMath::DegToRad();
    
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
    double m2 = M_eta;
    double m3 = M_He4;
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
    double P_LAB_eta = sqrt(pow(E_LAB_eta, 2) - pow(M_eta, 2));
    
    // Calculate the momentum for each direction
    double Px_LAB_eta = P_LAB_eta * sin(ThetaLAB) * cos(PhiLAB);
    double Py_LAB_eta = P_LAB_eta * sin(ThetaLAB) * sin(PhiLAB);
    double Pz_LAB_eta = P_LAB_eta * cos(ThetaLAB);
    
    // Store the results in TLorentzVector for the eta-meson
    TLorentzVector eta_LAB_4Vec(Px_LAB_eta, Py_LAB_eta, Pz_LAB_eta, E_LAB_eta);
    //TLorentzVector eta_COM_4Vec = eta_LAB_4Vec;
    //eta_COM_4Vec.Boost(-IS_4Vec.BoostVector());
        
    // Make the eta-meson decay into two photons
    double masses[] = {M_gamma, M_gamma};
    TLorentzVector photon_4Vec[2];
    if (decayGen.SetDecay(eta_LAB_4Vec, 2, masses)) {
      decayGen.Generate();
      photon_4Vec[0] = * decayGen.GetDecay(0);
      photon_4Vec[1] = * decayGen.GetDecay(1);
    }
    
    // Deduce by energy and mass conservation the recoil nucleus 4Vec
    TLorentzVector He4_LAB_4Vec = IS_4Vec - eta_LAB_4Vec;
    
    h_Tkin_recoilA_vs_egam->Fill(ebeam, He4_LAB_4Vec.E() - He4_LAB_4Vec.M());
    h_theta_recoilA_vs_egam->Fill(ebeam, He4_LAB_4Vec.Theta() * TMath::RadToDeg());
    h_Tkin_eta_vs_egam->Fill(ebeam, eta_LAB_4Vec.E() - eta_LAB_4Vec.M());
    h_theta_eta_vs_egam->Fill(ebeam, eta_LAB_4Vec.Theta() * TMath::RadToDeg());
    h_Tkin_photon_vs_egam->Fill(ebeam, photon_4Vec[0].E());
    h_theta_photon_vs_egam->Fill(ebeam, photon_4Vec[0].Theta() * TMath::RadToDeg());
    h_Tkin_photon_vs_egam->Fill(ebeam, photon_4Vec[1].E());
    h_theta_photon_vs_egam->Fill(ebeam, photon_4Vec[1].Theta() * TMath::RadToDeg());
    
    if (hddmWriter) {
      // ======= HDDM output =========
      tmpEvt_t tmpEvt;
      tmpEvt.beam = InGamma_4Vec;
      tmpEvt.target = Target_4Vec;
      tmpEvt.q1 = photon_4Vec[0];
      tmpEvt.q2 = photon_4Vec[1];
      tmpEvt.q3 = He4_LAB_4Vec;
      tmpEvt.nGen = 3;
      tmpEvt.weight = 1.;
      hddmWriter->write(tmpEvt,runNum,i);
    }
    if (asciiWriter) {
      // ======= ASCII output =========
      (*asciiWriter)<<runNum<<" "<<i<<" 3"<<endl;
      // photons from the eta
      (*asciiWriter)<<"0 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      (*asciiWriter)<<"   "<<0<<" "<<photon_4Vec[0].Px()<<" "<<photon_4Vec[0].Py()<<" "<<photon_4Vec[0].Pz()<<" "<<photon_4Vec[0].E()<<endl;			
      (*asciiWriter)<<"1 "<<gamma_TYPE<<" "<<M_gamma<<endl;
      (*asciiWriter)<<"   "<<0<<" "<<photon_4Vec[1].Px()<<" "<<photon_4Vec[1].Py()<<" "<<photon_4Vec[1].Pz()<<" "<<photon_4Vec[1].E()<<endl;			
      // Nucleus recoil
      (*asciiWriter)<<"2 "<<Helium_TYPE<<" "<<M_He4<<endl;
      (*asciiWriter)<<"   "<<1<<" "<<He4_LAB_4Vec.Px()<<" "<<He4_LAB_4Vec.Py()<<" "<<He4_LAB_4Vec.Pz()<<" "<<He4_LAB_4Vec.E()<<endl;			
    }
    
    // deletion
    delete h_ThetaLAB;
  }
  h_Tkin_eta_vs_egam->Write();
  h_Tkin_photon_vs_egam->Write();
  h_Tkin_recoilA_vs_egam->Write();
  h_theta_eta_vs_egam->Write();
  h_theta_photon_vs_egam->Write();
  h_theta_recoilA_vs_egam->Write();
  //h_dxs->Write();
  //cobrem_vs_Erange->Write();
  diagOut->Close();
  
  if (hddmWriter) delete hddmWriter;
  if (asciiWriter) delete asciiWriter;
  
  return 0;
}


