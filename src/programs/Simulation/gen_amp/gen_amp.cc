
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
#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"

#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/BreitWigner3body.h"
#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
#include "AMPTOOLS_AMPS/Lambda1520Angles.h"
#include "AMPTOOLS_AMPS/Lambda1520tdist.h"
#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/Ylm.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/IsobarAngles.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToNPartP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

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
	double lowMass = 0.2;
	double highMass = 3.0;

	double beamMaxE   = 12.0;
	double beamPeakE  = 9.0;
	double beamLowE   = 2.0;
	double beamHighE  = 12.0;

	int runNum = 30731;
	unsigned int seed = 0;

	double lowT = 0.0;
	double highT = 12.0;
	double slope = 6.0;

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
		if (arg == "-t"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  slope = atof( argv[++i] ); }
		if (arg == "-tmin"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  lowT = atof( argv[++i] ); }
		if (arg == "-tmax"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  highT = atof( argv[++i] ); }
		if (arg == "-d"){
			diag = true; }
		if (arg == "-f"){
			genFlat = true; }
		if (arg == "-h"){
			cout << endl << " Usage for: " << argv[0] << endl << endl;
			cout << "\t -c    <file>\t Config file" << endl;
			cout << "\t -o    <name>\t ROOT file output name" << endl;
			cout << "\t -hd   <name>\t HDDM file output name [optional]" << endl;
			cout << "\t -l    <value>\t Low edge of mass range (GeV) [optional]" << endl;
			cout << "\t -u    <value>\t Upper edge of mass range (GeV) [optional]" << endl;
			cout << "\t -n    <value>\t Minimum number of events to generate [optional]" << endl;
			cout << "\t -m  <value>\t Electron beam energy (or photon energy endpoint) [optional]" << endl;
                        cout << "\t -p  <value>\t Coherent peak photon energy [optional]" << endl;
                        cout << "\t -a  <value>\t Minimum photon energy to simulate events [optional]" << endl;
                        cout << "\t -b  <value>\t Maximum photon energy to simulate events [optional]" << endl;
			cout << "\t -r    <value>\t Run number assigned to generated events [optional]" << endl;
			cout << "\t -s    <value>\t Random number seed initialization [optional]" << endl;
			cout << "\t -t    <value>\t Momentum transfer slope [optional]" << endl;
			cout << "\t -tmin <value>\t Minimum momentum transfer [optional]" << endl;
			cout << "\t -tmax <value>\t Maximum momentum transfer [optional]" << endl;
			cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
			cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
			exit(1);
		}
	}
	
	if( configfile.size() == 0 || outname.size() == 0 ){
		cout << "No config file or output specificed:  run gen_amp -h for help" << endl;
		exit(1);
	}
	
	// open config file and be sure only one reaction is specified
	ConfigFileParser parser( configfile );
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	assert( cfgInfo->reactionList().size() == 1 );
	ReactionInfo* reaction = cfgInfo->reactionList()[0];
	
	// use particletype.h to convert reaction particle names
	vector<Particle_t> Particles;
	vector<double> childMasses;
	double threshold = 0;
	for (unsigned int i = 0; i < reaction->particleList().size(); i++){
	  Particle_t locEnum = ParticleEnum(reaction->particleList()[i].c_str());
	  // Beam particle is always photon
	  if (locEnum == 0 && i > 0)
	    cout << "ConfigFileParser WARNING:  unknown particle type \"" << reaction->particleList()[i] << "\"" << endl;
	  Particles.push_back(ParticleEnum(reaction->particleList()[i].c_str()));
	  if (i>1){
	    childMasses.push_back(ParticleMass(Particles[i]));
	    threshold += ParticleMass(Particles[i]);
	  }
	}

	// loop to look for resonance in config file
	// currently only one at a time is supported 
	const vector<ConfigFileLine> configFileLines = parser.getConfigFileLines();
	double resonance[]={1.0, 1.0};
	bool foundResonance = false;
	bool isKaonRecoil = false;
	bool isPionRecoil = false;
	for (vector<ConfigFileLine>::const_iterator it=configFileLines.begin(); it!=configFileLines.end(); it++) {
	  if ((*it).keyword() == "define") {
	      if ((*it).arguments()[0] == "rho" || (*it).arguments()[0] == "omega" || (*it).arguments()[0] == "phi" || (*it).arguments()[0] == "b1" || (*it).arguments()[0] == "a1" || (*it).arguments()[0] == "Lambda1520" || (*it).arguments()[0] == "X"){
	      if ( (*it).arguments().size() != 3 )
		continue;
	      resonance[0]=atof((*it).arguments()[1].c_str());
	      resonance[1]=atof((*it).arguments()[2].c_str());
	      cout << "Distribution seeded with resonance " << (*it).arguments()[0] << " : mass = " << resonance[0] << "GeV , width = " << resonance[1] << "GeV" << endl; 
	      if((*it).arguments()[0] == "Lambda1520")
		 isKaonRecoil = true;
	      foundResonance = true;
	      break;
	    }
	  }
	}
	if (!foundResonance)
	  cout << "ConfigFileParser WARNING:  no known resonance found, seed with mass = width = 1GeV" << endl; 

	// random number initialization (set to 0 by default)
	TRandom3* gRandom = new TRandom3();
	gRandom->SetSeed(seed);
	seed = gRandom->GetSeed();
	cout << "TRandom3 Seed : " << seed << endl;

	// setup AmpToolsInterface
	AmpToolsInterface::registerAmplitude( ThreePiAngles() );
	AmpToolsInterface::registerAmplitude( TwoPiAngles() );
	AmpToolsInterface::registerAmplitude( TwoPiAngles_amp() );
	AmpToolsInterface::registerAmplitude( TwoPSHelicity() );
	AmpToolsInterface::registerAmplitude( TwoPSAngles() );
	AmpToolsInterface::registerAmplitude( BreitWigner() );
	AmpToolsInterface::registerAmplitude( BreitWigner3body() );
	AmpToolsInterface::registerAmplitude( ThreePiAnglesSchilling() );
	AmpToolsInterface::registerAmplitude( Lambda1520Angles() );
	AmpToolsInterface::registerAmplitude( Lambda1520tdist() );
	AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
	AmpToolsInterface::registerAmplitude( Ylm() );
	AmpToolsInterface::registerAmplitude( Zlm() );
	AmpToolsInterface::registerAmplitude( IsobarAngles() );
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

	// generate over a range of mass
	// start with threshold or lowMass, whichever is higher
	GammaPToNPartP resProd;
	if(isKaonRecoil)
		resProd = GammaPToNPartP( threshold<lowMass ? lowMass : threshold, highMass, childMasses, ProductionMechanism::kKaon, type, slope, lowT, highT, seed, beamConfigFile );
	else if(isPionRecoil)
		resProd = GammaPToNPartP( threshold<lowMass ? lowMass : threshold, highMass, childMasses, ProductionMechanism::kPion, type, slope, lowT, highT, seed, beamConfigFile );
	else
		resProd = GammaPToNPartP( threshold<lowMass ? lowMass : threshold, highMass, childMasses, ProductionMechanism::kProton, type, slope, lowT, highT, seed, beamConfigFile );
	
	if (childMasses.size() < 2){
	  cout << "ConfigFileParser ERROR:  single particle production is not yet implemented" << endl; 
	  return 1;
	}
		
	// seed the distribution with a sum of noninterfering Breit-Wigners
	// we can easily compute the PDF for this and divide by that when
	// doing accept/reject -- improves efficiency if seeds are picked well
	
	if( !genFlat ){
		
		// the lines below should be tailored by the user for the particular desired
		// set of amplitudes -- doing so will improve efficiency.  Leaving as is
		// won't make MC incorrect, it just won't be as fast as it could be
		
		resProd.addResonance( resonance[0], resonance[1],  1.0 );
	}
	
	vector< int > pTypes;
	for (unsigned int i=0; i<Particles.size(); i++)
	  pTypes.push_back( Particles[i] );
	
	HDDMDataWriter* hddmOut = NULL;
	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, runNum, seed);
	ROOTDataWriter rootOut( outname );
	
	TFile* diagOut = new TFile( "gen_amp_diagnostic.root", "recreate" );
	ostringstream locStream;
	ostringstream locIsobarStream;
	for (unsigned int i=2; i<Particles.size(); i++){
	  locStream << ParticleName_ROOT(Particles[i]);
	  if ( i> 2 )
	    locIsobarStream << ParticleName_ROOT(Particles[i]);
	}
	string locHistTitle = string("Resonance Mass ;") + locStream.str() + string(" Invariant Mass (GeV/c^{2});");
	string locIsobarTitle = string("Isobar Mass ;") + locIsobarStream.str() + string(" Invariant Mass (GeV/c^{2});");

	TH1F* mass = new TH1F( "M", locHistTitle.c_str(), 180, lowMass, highMass );
	TH1F* massW = new TH1F( "M_W", ("Weighted "+locHistTitle).c_str(), 180, lowMass, highMass );
	massW->Sumw2();
	TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
	TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
	
	TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );

	TH1F* M_isobar = new TH1F( "M_isobar", locIsobarTitle.c_str(), 200, 0, 2 );

	TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
	TH2F* M_CosTheta = new TH2F( "M_CosTheta", "M vs. cos#vartheta", 180, lowMass, highMass, 200, -1, 1);
	TH2F* M_Phi = new TH2F( "M_Phi", "M vs. #varphi", 180, lowMass, highMass, 200, -3.14, 3.14);
	TH2F* M_Phi_lab = new TH2F( "M_Phi_lab", "M vs. #varphi", 180, lowMass, highMass, 200, -3.14, 3.14);

	TH2F* CosTheta_Phi = new TH2F( "CosTheta_Phi", "cos#theta vs. #phi", 180, -3.14, 3.14, 100, -1, 1);

	// define isobars and create histograms
	int nIsobars = 2;
	vector<int> isobarIndex[nIsobars];
	isobarIndex[0].push_back(2); 
	isobarIndex[0].push_back(3);
	//isobarIndex[0].push_back(4);
	isobarIndex[1].push_back(4);
	isobarIndex[1].push_back(5);
	
	TH2F* CosTheta_PhiIso[nIsobars];
	for(int i=0; i<nIsobars; i++)
		CosTheta_PhiIso[i] = new TH2F( Form("CosTheta_PhiIso%d",i), "cos#theta vs. #phi", 180, -3.14, 3.14, 100, -1, 1);

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
			TLorentzVector resonance;
			for (unsigned int i=2; i<Particles.size(); i++)
			  resonance += evt->particle( i );

			TLorentzVector isobar;
			for (unsigned int i=2; i<4 - 1; i++)
			  isobar += evt->particle( i );

			/////////////////////////////////////////////////////////////////////
			// Do user calculation here for decay angles to match IsobarAngles //
			/////////////////////////////////////////////////////////////////////
			TLorentzVector beam = evt->particle(0);
			TLorentzVector recoil = evt->particle(1);
			TLorentzVector PX = resonance;
			TLorentzVector PBatchX = resonance - isobar;
			
			// calculate decay angles in resonance X rest frame
			TVector3 XRestBoost = PX.BoostVector();
			
			TLorentzVector beamX   = beam;
			TLorentzVector recoilX = recoil;
			TLorentzVector batchX  = PBatchX;
			beamX.Boost(-1.0*XRestBoost);
			recoilX.Boost(-1.0*XRestBoost);
			batchX.Boost(-1.0*XRestBoost);
			
			TVector3 z = beamX.Vect().Unit(); //-recoilX.Vect().Unit();
			TVector3 y = (beamX.Vect().Unit()).Cross((-recoilX.Vect().Unit())).Unit();
			TVector3 x = y.Cross(z);
			
			TVector3 anglesBatchX( (batchX.Vect()).Dot(x),
					       (batchX.Vect()).Dot(y),
					       (batchX.Vect()).Dot(z) );
			
			GDouble cosThetaBatchX = anglesBatchX.CosTheta();
			GDouble phiBatchX = anglesBatchX.Phi();
			
			// build 4-vectors for isobars (only for subsequent isobars now...)
			TLorentzVector PIsobar[nIsobars], PBatch[nIsobars];
			TVector3 zIsoX[nIsobars];
			for(int i=0; i<nIsobars; i++) {
				cout<<"Isobar "<<i<<endl;
				for(uint j=0; j<isobarIndex[i].size(); j++) { 
					cout<<"  Particle "<<j<<" mass = "<<(evt->particle(isobarIndex[i][j])).M()<<endl;
					if(j==0) PBatch[i] = evt->particle(isobarIndex[i][j]); // this is the problem for non-subsequent bachelors..
					PIsobar[i] += evt->particle(isobarIndex[i][j]);
				}
				
				TLorentzVector PIsobarX = PIsobar[i];
				PIsobarX.Boost(-1.0*XRestBoost);
				zIsoX[i] = PIsobarX.Vect().Unit();
			}
				
			////////////////////////////////////////////////////////////////////////////
			// calculate decay angles in isobar rest frame (NEED TO CHECK FOR BUGS!)) //
			////////////////////////////////////////////////////////////////////////////
			vector<GDouble> cosThetaIso, phiIso;
			pair<TVector3, TVector3> zIsoPrevious;
			for( int i = 0; i < nIsobars; i++ ){
				
				TVector3 isoRestBoost = PIsobar[i].BoostVector();
				TLorentzVector PIsobarIso = PIsobar[i];
				TLorentzVector PBatchIso = PBatch[i];
				TLorentzVector PResonanceIso = PIsobar[i] - PBatch[i];
				//TLorentzVector PNormIso1 = PNorm[i].first;
				//TLorentzVector PNormIso2 = PNorm[i].second;
				PBatchIso.Boost(-1.0*isoRestBoost);
				PResonanceIso.Boost(-1.0*isoRestBoost);
				PIsobarIso.Boost(-1.0*isoRestBoost);
				//PNormIso1.Boost(-1.0*isoRestBoost);
				//PNormIso2.Boost(-1.0*isoRestBoost);
				
				TVector3 zIso = zIsoX[i]; // z-axis is direction of isobar in X	rest frame by default
				TVector3 yIso = z.Cross(zIso); // decay plane from X rest frame
				
				// later stage of single batchelor decays (eg. omega->3pi in b1pi production)
				if(0) {//i == 1) { //m_isBach[i]) {
					zIso = zIsoPrevious.first;
					yIso = zIsoPrevious.second.Cross(zIsoPrevious.first);
				}
				TVector3 xIso = yIso.Cross(zIso);
				
				TVector3 PAngles = PBatchIso.Vect();
				//if(m_daughtI[i].size() == 3 and m_daughtI.size() == uint(m_nIsobars)) // 3-body decay (e.g. omega) 
				//	PAngles = (PNormIso1.Vect()).Cross(PNormIso2.Vect());
				
				// Angles in isobar rest frame
				TVector3 anglesIso( (PAngles).Dot(xIso),
						    (PAngles).Dot(yIso),
						    (PAngles).Dot(zIso) );
				
				cosThetaIso.push_back(anglesIso.CosTheta());
				phiIso.push_back(anglesIso.Phi());
				cout<<i<<" "<<cosThetaIso[i]<<endl;

				// reference vector for later step in current frame
				zIsoPrevious.first = PAngles.Unit();
				// reference vector for later step in previous frame
				zIsoPrevious.second = zIso;
			}
			
			double genWeight = evt->weight();
			
			// cannot ask for the intensity if we haven't called process events above
			double weightedInten = ( genFlat ? 1 : ati.intensity( i ) ); 
			// cout << " i=" << i << "  intensity_i=" << weightedInten << endl;

			if( !diag ){
				
				// obtain this by looking at the maximum value of intensity * genWeight
				double rand = gRandom->Uniform() * maxInten;
				
				if( weightedInten > rand || genFlat ){

					mass->Fill( resonance.M() );
					massW->Fill( resonance.M(), genWeight );
					double loct = -1.*(TLorentzVector(0,0,0,0.938) - evt->particle(1)).M2();
					t->Fill(loct);

					intenW->Fill( weightedInten );
					intenWVsM->Fill( resonance.M(), weightedInten );

					M_isobar->Fill( isobar.M() );

					M_CosTheta->Fill( resonance.M(), cosThetaBatchX);
					M_Phi->Fill( resonance.M(), phiBatchX);
					M_Phi_lab->Fill( resonance.M(), recoil.Phi());
					
					CosTheta_Phi->Fill( phiBatchX, cosThetaBatchX);
					for(int i=0; i<nIsobars; i++) 
						CosTheta_PhiIso[i]->Fill( phiIso[i], cosThetaIso[i]);

					// we want to save events with weight 1
					evt->setWeight( 1.0 );
					
					if( hddmOut ) hddmOut->writeEvent( *evt, pTypes );
					rootOut.writeEvent( *evt );
					++eventCounter;
					if(eventCounter >= nEvents) break;
				}
			}
			else{
				
				mass->Fill( resonance.M() );
				massW->Fill( resonance.M(), genWeight );
				
				intenW->Fill( weightedInten );
				intenWVsM->Fill( resonance.M(), weightedInten );
				TLorentzVector recoil = evt->particle ( 1 );
				
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
	M_isobar->Write();
	t->Write();
	CosTheta_psi->Write();
	M_CosTheta->Write();
	M_Phi->Write();
	M_Phi_lab->Write();

	CosTheta_Phi->Write();
	for(int i=0; i<nIsobars; i++) CosTheta_PhiIso[i]->Write();

	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}

