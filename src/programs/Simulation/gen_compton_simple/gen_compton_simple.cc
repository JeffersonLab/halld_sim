
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

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

#include "HddmOut.h"

using std::complex;
using namespace std;

#define GAMMA_TYPE 1
#define ELECTRON_TYPE 3


double gen_costheta(double ebeam)
{
      const double M_e = 0.51099892e-3;
      double random_factor = 2.;
      double diff_xs = 1.;
      double costheta = -1.;
      
      while(random_factor > diff_xs) {
		  double e      = ebeam / M_e;
	      //costheta      = 2.*(gRandom->Uniform(-1.,1.)-0.5);
	      costheta      = gRandom->Uniform(-1.,1.);
		  double p      = 1./(1. + e * (1.-costheta));
		  diff_xs  = 0.5*p*(p*p + 1. - (1.-costheta)*(1.+costheta)*p);
		  random_factor = gRandom->Uniform(0.,1.);
	}
	
	return costheta;
}

void calculate_compton_kinematics(double ebeam, double gamma_costheta, double &gamma_en,
								  double &electron_costheta, double &electron_en)
{
      const double M_e = 0.51099892e-3;

      double e   = ebeam / M_e;
      gamma_en  = ebeam / (1. + e*(1.-gamma_costheta));
      electron_en  = ebeam + M_e - gamma_en;
      electron_costheta  = cos(atan(1./(1.+e)/tan(0.5*acos(gamma_costheta))));
}


int main( int argc, char* argv[] ){
  
	string  configfile("");
	string  outname("");
	string  hddmname("");
	
	double beamLowE   = 3.0;
	double beamHighE  = 12.0;	
	const double M_e = 0.51099892e-3;
	
	int runNum = 9001;
	int seed = 0;

	int nEvents = 10000;
		
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
	
	if( configfile.size() == 0 || outname.size() == 0 ){
	  cout << "No config file or output specificed:  run gen_compton -h for help" << endl;
	  exit(1);
	}
	
	if( outname.size() == 0 && hddmname == 0){
		cout << "No output specificed:  run gen_compton_simple -h for help" << endl;
		exit(1);
	}

	// random number initialization (set to 0 by default)
	gRandom->SetSeed(seed);
	
	if( outname.size() == 0 && hddmname == 0){
		cout << "No output specificed:  run gen_compton_simple -h for help" << endl;
		exit(1);
	}
	
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
	TString beamConfigFile;
	BeamProperties beamProp(beamConfigFile);
	TH1D * cobrem_vs_E = (TH1D*)beamProp.GetFlux();
		
	for( int i = 0; i < nEvents; ++i ) {
		if(i%1000 == 1)
			cout << "event " << i <<endl;
		
		// get beam energy
		double ebeam = 0;
		  if (beamConfigFile != "")
		  ebeam = cobrem_vs_E->GetRandom();
		if (beamConfigFile == "")
		  ebeam = ebeam_spectrum.GetRandom();
		// generate cos(theta) according to dsig/dOmega
		double gamma_costheta = gen_costheta(ebeam);

		double electron_en = 0.;
		double gamma_en = 0.;
		double electron_costheta = 0.;
		
		// calculate the rest of the kinematics
		calculate_compton_kinematics(ebeam, gamma_costheta, gamma_en,
									 electron_costheta, electron_en);
									 
		// more kinematics
		double electron_momentum = sqrt(electron_en*electron_en - M_e*M_e);
		double gamma_theta = acos(gamma_costheta);
		double electron_theta = acos(electron_costheta);

		double electron_phi = gRandom->Uniform(0., 2.*TMath::Pi());
		double gamma_phi = electron_phi - TMath::Pi();
		if(gamma_phi < 0.)
			gamma_phi += 2.*TMath::Pi();

		TLorentzVector electron_4v(electron_momentum*cos(electron_phi)*sin(electron_theta), 
							  	   electron_momentum*sin(electron_phi)*sin(electron_theta), 
							  	   electron_momentum*electron_costheta, electron_en);
		TLorentzVector gamma_4v(gamma_en*cos(gamma_phi)*sin(gamma_theta), 
							    gamma_en*sin(gamma_phi)*sin(gamma_theta), 
							    gamma_en*gamma_costheta, gamma_en);
		TLorentzVector beam_4v(0., 0., ebeam, ebeam);
		TLorentzVector target_4v(0., 0., 0., 0.938);
		
		
		if(hddmWriter) {
			// ======= HDDM output =========
			tmpEvt_t tmpEvt;
			tmpEvt.beam = beam_4v;
			tmpEvt.target = target_4v;
			tmpEvt.q1 = electron_4v;
			tmpEvt.q2 = gamma_4v;
			//tmpEvt.recoil = proton_4v;
			tmpEvt.nGen = 2;
			tmpEvt.weight = 1.;
			hddmWriter->write(tmpEvt,runNum,i);
		}
		if(asciiWriter) {
			// ======= ASCII output =========
			(*asciiWriter)<<runNum<<" "<<i<<" 2"<<endl;
			// compton photon
			(*asciiWriter)<<"0 "<<GAMMA_TYPE<<" "<<0<<endl;
			(*asciiWriter)<<"   "<<0<<" "<<gamma_4v.Px()<<" "<<gamma_4v.Py()<<" "<<gamma_4v.Pz()<<" "<<gamma_4v.E()<<endl;			
			// compton electron
			(*asciiWriter)<<"1 "<<ELECTRON_TYPE<<" "<<0.000511<<endl;
			(*asciiWriter)<<"   "<<-1<<" "<<electron_4v.Px()<<" "<<electron_4v.Py()<<" "<<electron_4v.Pz()<<" "<<electron_4v.E()<<endl;;			
		}
	}
	
	
	if( hddmWriter ) delete hddmWriter;
	if( asciiWriter ) delete asciiWriter;
	
	return 0;
}


