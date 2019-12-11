
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

#include "HddmOut.h"

using std::complex;
using namespace std;

#define eta_TYPE 17
#define Helium4_TYPE 1000020040

int main( int argc, char* argv[] ){
  
	string  configfile("");
	string  outname("");
	string  hddmname("");
	
	ifstream in_coherent;
	ifstream in_incoherent;

	double beamLowE   = 3.0;
	double beamHighE  = 12.0;
	const double M_He4 = 4.002602;
	const double M_eta = 0.54730;
	TLorentzVector Target_4Vec(0, 0, 0, M_He4);
	
	int runNum = 9001;
	int seed = 0;

	int nEvents = 10000;
	
	double dummy;
	
	//parse command line:
	for (int i = 1; i < argc; i++){
		
		string arg(argv[i]);
		
		if (arg == "-c"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  configfile = argv[++i]; }
		if (arg == "-o"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  outname = argv[++i]; }
		if (arg == "-hd"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  hddmname = argv[++i]; }
		if (arg == "-n"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  nEvents = atoi( argv[++i] ); }
		if (arg == "-a"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamLowE = atof( argv[++i] ); }
		if (arg == "-b"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamHighE = atof( argv[++i] ); }
		if (arg == "-r"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  runNum = atoi( argv[++i] ); }
		if (arg == "-s"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  seed = atoi( argv[++i] ); }
		if (arg == "-h"){
			cout << endl << " Usage for: " << argv[0] << endl << endl;
			cout << "\t -c  <file>\t Config file" << endl;
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
	
	if( outname.size() == 0 && hddmname == 0){
		cout << "No output specificed:  run gen_compton_simple -h for help" << endl;
		exit(1);
	}
	
	// random number initialization (set to 0 by default)
	gRandom->SetSeed(seed);

	// initialize HDDM output
	HddmOut *hddmWriter = nullptr;
	if(hddmname != "")
		hddmWriter = new HddmOut(hddmname.c_str());
		
	// initialize ASCII output
	ofstream *asciiWriter = nullptr;
	if(outname != "")
		asciiWriter = new ofstream(outname.c_str());

	// Assume a beam energy spectrum of 1/E(gamma)
	TF1 ebeam_spectrum("beam_spectrum","1/x",beamLowE,beamHighE);
	
	// get beam properties from configuration file
	TH1D * cobrem_vs_E = 0;
	if (configfile != "") {
	  BeamProperties beamProp( configfile );
	  cobrem_vs_E = (TH1D*)beamProp.GetFlux();
	}
	for (int i = 0; i < nEvents; ++i) {
		if(i%1000 == 1)
		  cout << "event " << i <<endl;
	
		// get beam energy
		double ebeam = 0;
		if (configfile == "") 
		  ebeam = ebeam_spectrum.GetRandom();
		else if (configfile != "")
		  ebeam = cobrem_vs_E->GetRandom();
		
		// Incident photon-beam 4Vec
		TLorentzVector InGamma_4Vec(0, 0, ebeam, ebeam);
		
		// Initial state 4Vec
		TLorentzVector IS_4Vec = InGamma_4Vec + Target_4Vec;
		
		// Mass in the centre-of-mass frame
		double sqrt_s = IS_4Vec.M();
		double s = pow(sqrt_s, 2);
		
		// Histo. creation that will store the calculated diff. xs. vs. LAB polar angle
		TH1F * h_ThetaLAB = new TH1F("h_ThetaLAB", "", 300, 0., 3.0);
		
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
		
		// Generate eta-meson theta in LAB
		double ThetaLAB = h_ThetaLAB->GetRandom();

		// Generate eta-meson phi in LAB
		double PhiLAB = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
		
		// Calculate eta-meson energy in COM
		double E_COM_eta = (s - pow(M_He4, 2) + pow(M_eta, 2)) / (2.0 * sqrt_s);
		
		// Calculate eta-meson momentum in COM
		double P_COM_eta = sqrt(pow(E_COM_eta, 2) - pow(M_eta, 2));

		// Calculate eta-meson momentun im LAB
		double P_LAB_eta = P_COM_eta * sqrt_s / M_He4;
		
		// Calculate eta-meson energy in LAB
		double E_LAB_eta = sqrt(pow(P_LAB_eta, 2) + pow(M_eta, 2));
		
		// Calculate the momentum for each direction
		double Px_LAB_eta = P_LAB_eta * sin(ThetaLAB) * cos(PhiLAB);
		double Py_LAB_eta = P_LAB_eta * sin(ThetaLAB) * sin(PhiLAB);
		double Pz_LAB_eta = P_LAB_eta * cos(ThetaLAB);
		
		// Store the results in TLorentzVector for the eta-meson
		TLorentzVector eta_LAB_4Vec(Px_LAB_eta, Py_LAB_eta, Pz_LAB_eta, E_LAB_eta);
		
		// Deduce by energy and mass conservation the recoil nucleus 4Vec
		TLorentzVector He4_LAB_4Vec = IS_4Vec - eta_LAB_4Vec;
		
		if(hddmWriter) {
			// ======= HDDM output =========
			tmpEvt_t tmpEvt;
			tmpEvt.beam = InGamma_4Vec;
			tmpEvt.target = Target_4Vec;
			tmpEvt.q1 = eta_LAB_4Vec;
			tmpEvt.q2 = He4_LAB_4Vec;
			tmpEvt.nGen = 2;
			tmpEvt.weight = 1.;
			hddmWriter->write(tmpEvt,runNum,i);
		}
		if(asciiWriter) {
			// ======= ASCII output =========
			(*asciiWriter)<<runNum<<" "<<i<<" 2"<<endl;
			// compton photon
			(*asciiWriter)<<"0 "<<eta_TYPE<<" "<<M_eta<<endl;
			(*asciiWriter)<<"   "<<0<<" "<<eta_LAB_4Vec.Px()<<" "<<eta_LAB_4Vec.Py()<<" "<<eta_LAB_4Vec.Pz()<<" "<<eta_LAB_4Vec.E()<<endl;			
			// compton electron
			(*asciiWriter)<<"1 "<<Helium4_TYPE<<" "<<M_He4<<endl;
			(*asciiWriter)<<"   "<<-1<<" "<<He4_LAB_4Vec.Px()<<" "<<He4_LAB_4Vec.Py()<<" "<<He4_LAB_4Vec.Pz()<<" "<<He4_LAB_4Vec.E()<<endl;;			
		
		}
		
		// deletion
		delete h_ThetaLAB;
	}
	
	
	if( hddmWriter ) delete hddmWriter;
	if( asciiWriter ) delete asciiWriter;
	
	return 0;
}


