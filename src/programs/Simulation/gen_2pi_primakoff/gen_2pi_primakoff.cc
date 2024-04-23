
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

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "AMPTOOLS_MCGEN/HDDMDataWriter.h"

#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/TwoPiWt_primakoff.h"
#include "AMPTOOLS_AMPS/TwoPiWt_sigma.h"
#include "AMPTOOLS_AMPS/TwoPitdist.h"
#include "AMPTOOLS_AMPS/TwoPiNC_tdist.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaZToXYZ.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
  
	string  configfile("");
	string  outname("");
	string  hddmname("");
	
	bool diag = false;
	bool genFlat = false;
	
	// default upper and lower bounds 
	// double lowMass = 0.2;
	// double highMass = 2.0; 
	double PiMass = ParticleMass(PiPlus);
	double lowMass = 2*PiMass + 0.01;   // slightly higher than actual threshold to avoid zeros.;              
	double highMass = 0.8 ;
	
	double beamMaxE   = 12.0;
	double beamPeakE  = 6.0;
	double beamLowE   = lowMass;
	double beamHighE  = 12.0;
	
	int runNum = 9001;
	int seed = 0;

	int nEvents = 10000;
	int batchSize = 10000;
	
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
		if (arg == "-l"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  lowMass = atof( argv[++i] ); }
		if (arg == "-u"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  highMass = atof( argv[++i] ); }
		if (arg == "-n"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  nEvents = atoi( argv[++i] ); }
		if (arg == "-m"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamMaxE = atof( argv[++i] ); }
		if (arg == "-p"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else beamPeakE = atof( argv[++i] ); }
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
		if (arg == "-d"){
			diag = true; }
		if (arg == "-f"){
			genFlat = true; }
		if (arg == "-h"){
			cout << endl << " Usage for: " << argv[0] << endl << endl;
			cout << "\t -c  <file>\t Config file" << endl;
			cout << "\t -o  <name>\t ROOT file output name" << endl;
			cout << "\t -hd <name>\t HDDM file output name [optional]" << endl;
			cout << "\t -l  <value>\t Low edge of mass range (GeV) [optional]" << endl;
			cout << "\t -u  <value>\t Upper edge of mass range (GeV) [optional]" << endl;
			cout << "\t -n  <value>\t Minimum number of events to generate [optional]" << endl;
			cout << "\t -m  <value>\t Electron beam energy (or photon energy endpoint) [optional]" << endl;
                        cout << "\t -p  <value>\t Coherent peak photon energy [optional]" << endl;
                        cout << "\t -a  <value>\t Minimum photon energy to simulate events [optional]" << endl;
                        cout << "\t -b  <value>\t Maximum photon energy to simulate events [optional]" << endl;
			cout << "\t -r  <value>\t Run number assigned to generated events [optional]" << endl;
			cout << "\t -s  <value>\t Random number seed initialization [optional]" << endl;
			cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
			cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
			exit(1);
		}
	}
	
	if( configfile.size() == 0 || outname.size() == 0 ){
		cout << "No config file or output specificed:  run gen_2pi_primakoff -h for help" << endl;
		exit(1);
	}
	
	// open config file and be sure only one reaction is specified
	ConfigFileParser parser( configfile );
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	assert( cfgInfo->reactionList().size() == 1 );
	ReactionInfo* reaction = cfgInfo->reactionList()[0];
	
	// random number initialization (set to 0 by default)
	gRandom->SetSeed(seed);

	// setup AmpToolsInterface
	AmpToolsInterface::registerAmplitude( TwoPiAngles_primakoff() );
	AmpToolsInterface::registerAmplitude( TwoPiWt_primakoff() );
	AmpToolsInterface::registerAmplitude( TwoPiWt_sigma() );
	AmpToolsInterface::registerAmplitude( TwoPitdist() );	
	AmpToolsInterface::registerAmplitude( TwoPiNC_tdist() );
	AmpToolsInterface::registerAmplitude( BreitWigner() );
	AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );

	// loop to look for beam configuration file
        TString beamConfigFile;
        const vector<ConfigFileLine> configFileLinesBeam = parser.getConfigFileLines();
        for (vector<ConfigFileLine>::const_iterator it=configFileLinesBeam.begin(); it!=configFileLinesBeam.end(); it++) {
                if ((*it).keyword() == "define") {
                        TString beamArgument =  (*it).arguments()[0].c_str();
                        if(beamArgument.Contains("beamconfig")) {
                                beamConfigFile = (*it).arguments()[1].c_str();
                        }
                }
        }
	if(beamConfigFile.Length() == 0) {
		cout<<"WARNING: Couldn't find beam configuration file -- write local version"<<endl;

		beamConfigFile = "local_beam.conf";
		ofstream locBeamConfigFile;
		locBeamConfigFile.open(beamConfigFile.Data());
		locBeamConfigFile<<"ElectronBeamEnergy "<<beamMaxE<<endl;       // electron beam energy
		locBeamConfigFile<<"CoherentPeakEnergy "<<beamPeakE<<endl;      // coherent peak energy
		locBeamConfigFile<<"PhotonBeamLowEnergy "<<beamLowE<<endl;      // photon beam low energy
		locBeamConfigFile<<"PhotonBeamHighEnergy "<<beamHighE<<endl;    // photon beam high energy
		locBeamConfigFile.close();
	}
	
	ProductionMechanism::Type type =
		( genFlat ? ProductionMechanism::kFlat : ProductionMechanism::kResonant );
	
	// generate over a range of mass -- the daughters are two charged pions
	float Bgen= 230;   // exponential slope, make it smaller than any slope in the generator.
	GammaZToXYZ resProd( lowMass, highMass, PiMass, PiMass, type, beamConfigFile , Bgen);
	
	// seed the distribution with a sum of noninterfering s-wave amplitudes
	// we can easily compute the PDF for this and divide by that when
	// doing accept/reject -- improves efficiency if seeds are picked well
	
	if( !genFlat ){
		
		// the lines below should be tailored by the user for the particular desired
		// set of amplitudes -- doing so will improve efficiency.  Leaving as is
		// won't make MC incorrect, it just won't be as fast as it could be
		
		// resProd.addResonance( 0.775, 0.146,  1.0 );
		resProd.addResonance( 0.4, 0.146,  1.0 );    
	}
	
	vector< int > pTypes;
	pTypes.push_back( Gamma );
	pTypes.push_back( PiPlus );
	pTypes.push_back( PiMinus );
	pTypes.push_back( Pb208 );     // use lead instead of Sn116 since it is defined in particle list.
	
	HDDMDataWriter* hddmOut = NULL;
	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, runNum );
	ROOTDataWriter rootOut( outname );
	
	TFile* diagOut = new TFile( "gen_2pi_primakoff_diagnostic.root", "recreate" );
	
	TH1F* mass = new TH1F( "M", "Resonance Mass", 180, lowMass, highMass );
	TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, lowMass, highMass );
	massW->Sumw2();
	TH1D* intenW = new TH1D( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
	intenW->SetCanExtend(TH1::kXaxis);
	TH2D* intenWVsM = new TH2D( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
	intenWVsM->SetCanExtend(TH2::kYaxis);
	TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
	
	TH1D* h1_phi = new TH1D( "h1_phi", "#phi", 180, -PI,PI );
	TH1D* h1_psi = new TH1D( "h1_psi", "#psi", 180, -PI,PI );
	TH1D* h1_Phi = new TH1D( "h1_Phi", "#Phi", 180, -PI,PI );

	int eventCounter = 0;
	while( eventCounter < nEvents ){
		
		if( batchSize < 1E4 ){
			
			cout << "WARNING:  small batches could have batch-to-batch variations\n"
			     << "          due to different maximum intensities!" << endl;
		}
		
		cout << "Generating four-vectors..." << endl;
		
		ati.clearEvents();
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* kin = resProd.generate();
			ati.loadEvent( kin, i, batchSize );
			delete kin;
		}
		
		cout << "Processing events..." << endl;
		
		// include factor of 1.5 to be safe in case we miss peak -- avoid
		// intensity calculation of we are generating flat data
		double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
		
		
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* evt = ati.kinematics( i );
			TLorentzVector resonance( evt->particle( 1 ) + 
                                  evt->particle( 2 ) );
			
			double genWeight = evt->weight();
			
			// cannot ask for the intensity if we haven't called process events above
			// double ResM = resonance.M();
		        // double intensity_i = ati.intensity( i );
			double weightedInten = ( genFlat ? 1 : ati.intensity( i ) ); 
			// cout << " i=" << i << "  intensity_i=" << intensity_i << " maxInten=" << maxInten << " ResM=" << ResM << endl;

			if( !diag ){
				
				// obtain this by looking at the maximum value of intensity * genWeight
				double rand = gRandom->Uniform() * maxInten;
				
				if( weightedInten > rand || genFlat ){
					
					mass->Fill( resonance.M() );
					massW->Fill( resonance.M(), genWeight );
					
					intenW->Fill( weightedInten );
					intenWVsM->Fill( resonance.M(), weightedInten );
					
					// calculate angular variables
					TLorentzVector beam = evt->particle ( 0 );
					TLorentzVector recoil = evt->particle ( 3 );
					TLorentzVector p1 = evt->particle ( 1 );
					TLorentzVector p2 = evt->particle ( 2 );

					if (isfinite(recoil.M()) && isfinite(p1.M())) {
					  // check for nan values in vectors

					  // cout << endl << " gen_2pi_primakoff particles " << " Mbeam=" << beam.M() << " Mrecoil=" << recoil.M() << " Mp1=" << p1.M() << endl;
					  // beam.Print(); recoil.Print(); p1.Print(); p2.Print(); resonance.Print();
			     
					Double_t phipol=0;    // hardwire angle of photon polarization in lab.
					TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab
					TLorentzRotation resonanceBoost( -resonance.BoostVector() );
					
					TLorentzVector beam_res = resonanceBoost * beam;
					TLorentzVector recoil_res = resonanceBoost * recoil;
					TLorentzVector p1_res = resonanceBoost * p1;

					// choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is defined as missing P4
					TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();   
					TVector3 z = -1. * recoil_res.Vect().Unit();
					TVector3 x = y.Cross(z).Unit();
					TVector3 angles( (p1_res.Vect()).Dot(x),
							 (p1_res.Vect()).Dot(y),
							 (p1_res.Vect()).Dot(z) );

					GDouble CosTheta = angles.CosTheta();
					GDouble phi = angles.Phi();
					// GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
					// GDouble sin2Theta = sin(2.*angles.Theta());

					GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

					GDouble psi = Phi - phi;               // define angle difference 
					if(psi < -1*PI) psi += 2*PI;
					if  (psi > PI) psi -= 2*PI;

                                        // double phi = angles.Phi();

					/*cout << endl << " gen_2pi_primakoff " << endl;
                                        cout << " Phi=" << Phi << endl;
					cout << " phi= " << phi << endl;
					cout << " psi=" << psi << endl;*/

					h1_phi->Fill(phi);
					h1_psi->Fill(psi);
					h1_Phi->Fill(Phi);
					CosTheta_psi->Fill( psi, CosTheta);
					
					// we want to save events with weight 1
					evt->setWeight( 1.0 );
					float vx = 0;
					float vy = 0;
					float vz = 1;   // vertex for CCP experiment
					
					if( hddmOut ) hddmOut->writeEvent( *evt, pTypes, vx, vy, vz);
					// note that there is no provision currently for vertex output in root file
					rootOut.writeEvent( *evt );
					++eventCounter;
					}
				}
			}
			else{
				
				mass->Fill( resonance.M() );
				massW->Fill( resonance.M(), genWeight );
				
				intenW->Fill( weightedInten );
				intenWVsM->Fill( resonance.M(), weightedInten );
				TLorentzVector recoil = evt->particle ( 3 );
				
				++eventCounter;
			}
			
			delete evt;
		}
		
		cout << eventCounter << " events were processed." << endl;
	}
	
	mass->Write();
	massW->Write();
	intenW->Write();
	intenWVsM->Write();
	CosTheta_psi->Write();
	h1_phi->Write();
	h1_psi->Write();
	h1_Phi->Write();
	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}


