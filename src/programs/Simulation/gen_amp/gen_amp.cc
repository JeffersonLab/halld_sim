
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
#include "AMPTOOLS_AMPS/Uniform.h"

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
	double highMass = 5.0;

	double beamMaxE   = 12.0;
	double beamPeakE  = 9.0;
	double beamLowE   = 3.0;
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
	    //cout << "Particle " << i << " Mass " << Particles[i]; 
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
	    if ((*it).arguments()[0] == "X" || (*it).arguments()[0] == "rho" || (*it).arguments()[0] == "omega" || (*it).arguments()[0] == "phi" || (*it).arguments()[0] == "b1" || (*it).arguments()[0] == "a1" || (*it).arguments()[0] == "Lambda1520"){
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
	AmpToolsInterface::registerAmplitude( Uniform() );
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

	double targetMass = ParticleMass(ParticleEnum("Proton"));
	if(recoil == ProductionMechanism::kZ)
		targetMass = ParticleMass(Particles[1]);
	double recMass = ParticleMass(Particles[1]);
	double cmEnergy = sqrt(targetMass*(targetMass + 2*beamLowE));
	if ( cmEnergy < minMass + recMass ){
	  cout << "ConfigFileParser ERROR:  Minimum photon energy not high enough to create resonance!" << endl;
	  return 1;
	}
	else if ( cmEnergy < highMass + recMass )
	  cout << "ConfigFileParser WARNING:  Minimum photon energy not high enough to guarantee flat mass distribution!" << endl;	
		
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
	for (unsigned int i=0; i<ParticlesLowerVertex.size(); i++) {
	  if(ParticlesLowerVertex[i] == Proton || ParticlesLowerVertex[i] == Neutron) continue;
          pTypes.push_back( ParticlesLowerVertex[i] );
	}

	HDDMDataWriter* hddmOut = NULL;
	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, runNum, seed);
	ROOTDataWriter rootOut( outname );
	
	TFile* diagOut = new TFile( "gen_amp_diagnostic.root", "recreate" );
	ostringstream locStream;
	ostringstream locIsobarStream;
	for (unsigned int i=2; i<Particles.size(); i++){
	  locStream << ParticleName_ROOT(Particles[i]);
	  if ( i> 3 )
	    locIsobarStream << ParticleName_ROOT(Particles[i]);
	}
	string locHistTitle = string("Resonance Mass ;") + locStream.str() + string(" Invariant Mass (GeV/c^{2});");
	string locIsobarTitle = string("Isobar Mass ;") + locIsobarStream.str() + string(" Invariant Mass (GeV/c^{2});");

	TH1F* mass = new TH1F( "M", locHistTitle.c_str(), 400, 0.0, highMass );
	TH1F* massW = new TH1F( "M_W", ("Weighted "+locHistTitle).c_str(), 180, 0.0, highMass );
	massW->Sumw2();
	TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
	TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, 0.0, highMass, 1000, 0, 10 );
	
	TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );

	TH1F* M_isobar = new TH1F( "M_isobar", locIsobarTitle.c_str(), 200, 0, 2 );

	TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
	TH2F* M_CosTheta = new TH2F( "M_CosTheta", "M vs. cos#vartheta", 400, 0.0, highMass, 200, -1, 1);
	TH2F* M_Phi = new TH2F( "M_Phi", "M vs. #varphi", 400, 0.0, highMass, 200, -3.14, 3.14);
	TH2F* M_Phi_lab = new TH2F( "M_Phi_lab", "M vs. #varphi", 400, 0.0, highMass, 200, -3.14, 3.14);

	TH2F* M_CosTheta_Iso[4];
	TH2F* M_Phi_Iso[4];
	TH2F* M_CosTheta_Iso1sq[4][3];
        TH2F* M_Phi_Iso1sq[4][3];
	TH2F* M_CosTheta_Iso2sq[4][3][2];
        TH2F* M_Phi_Iso2sq[4][3][2];
	for(int i=0; i<4; i++) {
		M_CosTheta_Iso[i] = new TH2F( Form("M_CosTheta_Iso_%d", i), "M vs. cos#vartheta", 400, 0.0, highMass, 200, -1, 1);
		M_Phi_Iso[i] = new TH2F( Form("M_Phi_Iso_%d", i), "M vs. #varphi", 400, 0.0, highMass, 200, -3.14, 3.14);

        	for(int j=0; j<3; j++) {
                	M_CosTheta_Iso1sq[i][j] = new TH2F( Form("M_CosTheta_Iso1_%d_%d", i, j), "M vs. cos#vartheta", 400, 0.0, highMass, 200, -1, 1);
                	M_Phi_Iso1sq[i][j] = new TH2F( Form("M_Phi_Iso1_%d_%d", i, j), "M vs. #varphi", 400, 0.0, highMass, 200, -3.14, 3.14);

			for(int k=0; k<2; k++) {
                        	M_CosTheta_Iso2sq[i][j][k] = new TH2F( Form("M_CosTheta_Iso2_%d_%d_%d", i, j, k), "M vs. cos#vartheta", 400, 0.0, highMass, 200, -1, 1);
                        	M_Phi_Iso2sq[i][j][k] = new TH2F( Form("M_Phi_Iso2_%d_%d_%d", i, j, k), "M vs. #varphi", 400, 0.0, highMass, 200, -3.14, 3.14);
			}
                }
	}

	//TH2F* CosTheta_psi_Iso1 = new TH2F( "CosTheta_psi_Iso1", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
	TH2F* M_CosTheta_Iso1 = new TH2F( "M_CosTheta_Iso1", "M vs. cos#vartheta _1", 400, 0.0, highMass, 200, -1, 1);
	TH2F* M_Phi_Iso1 = new TH2F( "M_Phi_Iso1", "M vs. #varphi _1", 400, 0.0, highMass, 200, -3.14, 3.14);
	TH2F* M_Phi_lab_Iso1 = new TH2F( "M_Phi_lab_Iso1", "M vs. #varphi _1lab", 400, 0.0, highMass, 200, -3.14, 3.14);
	TH2F* BeamE_CosTheta_Iso1 = new  TH2F( "BeamE_CosTheta_Iso1", "BeamE vs. cos#vartheta _1", 400, beamLowE, beamHighE, 200, -1, 1);
	TH2F* M_CosTheta_Diff_Iso1 = new TH2F( "M_CosTheta_Diff_Iso1", "M vs. cos#vartheta - cos#vartheta _1", 400, 0.0, highMass, 200, -2, 2);
	TH2F* M_Phi_Diff_Iso1 = new TH2F( "M_Phi_Diff_Iso1", "M vs. #varphi - #varphi _1", 400, 0.0, highMass, 200, -3.14 * 2, 3.14 * 2);

	TH2F* M_CosTheta_Iso2 = new TH2F( "M_CosTheta_Iso2", "M vs. cos#vartheta _2", 400, 0.0, highMass, 200, -1, 1);
	TH2F* M_Phi_Iso2 = new TH2F( "M_Phi_Iso2", "M vs. #varphi _2", 400, 0.0, highMass, 200, -3.14, 3.14);
	TH2F* M_Phi_lab_Iso2 = new TH2F( "M_Phi_lab_Iso2", "M vs. #varphi _2lab", 400, 0.0, highMass, 200, -3.14, 3.14);
	TH2F* M_CosTheta_Diff_Iso2 = new TH2F( "M_CosTheta_Diff_Iso2", "M vs. cos#vartheta - cos#vartheta _2", 400, 0.0, highMass, 200, -2, 2);
	TH2F* M_Phi_Diff_Iso2 = new TH2F( "M_Phi_Diff_Iso2", "M vs. #varphi - #varphi _2", 400, 0.0, highMass, 200, -3.14 * 2, 3.14 * 2);

	TH2F* ThetaCorr_Iso1 = new TH2F("ThetaCorr_Iso1", "; #theta 1; #theta 1INV", 400, -3.14, 3.14, 400, -3.14, 3.14);
	TH2F* PhiCorr_Iso1 = new TH2F("PhiCorr_Iso1", "; #phi 1; #phi 1INV", 400, -3.14*2, 3.14*2, 400, -3.14*2, 3.14*2);
	TH2F* ThetaCorr_Iso2 = new TH2F("ThetaCorr_Iso2", "; #theta 2; #theta 2INV", 400, -3.14, 3.14, 400, -3.14, 3.14);
	TH2F* PhiCorr_Iso2 = new TH2F("PhiCorr_Iso2", "; #phi 2; #phi 2INV", 400, -3.14*2, 3.14*2, 400, -3.14*2, 3.14*2);
	TH2F* PVsTheta[6];
	for(int ipart=1; ipart<6; ipart++) 
		PVsTheta[ipart] = new TH2F(Form("PVsTheta_%d",ipart),Form("Particle %d; #theta; p",ipart),100,0,3.14,120,0,12);

	int eventCounter = 0;
	while( eventCounter < nEvents ){
		
		if( batchSize < 1E4 ){
			
			cout << "WARNING:  small batches could have batch-to-batch variations\n"
			     << "          due to different maximum intensities!" << endl;
		}
		
		cout << "Generating four-vectors..." << endl;
		
		ati.clearEvents();
		int i=0;
                while( i < batchSize ){

			Kinematics* kin;
			if(bwGenLowerVertex.size() == 0) 
				kin = resProd.generate(); // stable particle at lower vertex
			else { 
				// unstable particle at lower vertex
				pair< double, double > bwLowerVertex = bwGenLowerVertex[0]();
				double lowerVertex_mass_bw = bwLowerVertex.first;
				if ( lowerVertex_mass_bw < thresholdLowerVertex || lowerVertex_mass_bw > 2.0) continue;
				resProd.getProductionMechanism().setRecoilMass( lowerVertex_mass_bw );
				
				Kinematics* step1 = resProd.generate();
				TLorentzVector beam = step1->particle( 0 );
				TLorentzVector recoil = step1->particle( 1 );
				
				// loop over meson decay
				vector<TLorentzVector> mesonChild;
				for(unsigned int i=0; i<childMasses.size(); i++) 
					mesonChild.push_back(step1->particle( 2+i ));
				
				// decay step for lower vertex
				TLorentzVector nucleon; // proton or neutron
				NBodyPhaseSpaceFactory lowerVertex_decay = NBodyPhaseSpaceFactory( lowerVertex_mass_bw, massesLowerVertex);
				vector<TLorentzVector> lowerVertexChild = lowerVertex_decay.generateDecay();
				// boost to lab frame via recoil kinematics
				for(unsigned int j=0; j<lowerVertexChild.size(); j++) 
				  lowerVertexChild[j].Boost( recoil.BoostVector() );
				nucleon = lowerVertexChild[0];

				// store particles in kinematic class
				vector< TLorentzVector > allPart;
				allPart.push_back( beam );
				allPart.push_back( nucleon );
				// loop over meson decay particles
				for(unsigned int j=0; j<mesonChild.size(); j++) 
					allPart.push_back(mesonChild[j]);
				// loop over lower vertex decay particles
				for(unsigned int j=1; j<lowerVertexChild.size(); j++) 
					allPart.push_back(lowerVertexChild[j]);
				
				kin = new Kinematics( allPart, 1.0 );
				delete step1;				
			}
			
			Kinematics* kin = resProd.generate();
			ati.loadEvent( kin, i, batchSize );
			delete kin;
			i++;
		}
		
		cout << "Processing events..." << endl;
		
		// include factor of 1.5 to be safe in case we miss peak -- avoid
		// intensity calculation of we are generating flat data
		//cout << "Particles.size() = " << Particles.size << endl; 
		double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
		
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* evt = ati.kinematics( i );
			TLorentzVector resonance;
			for (unsigned int i=2; i<Particles.size(); i++) 
			  resonance += evt->particle( i );

			TLorentzVector isobar;
			//for (unsigned int i=4; i<Particles.size(); i++)
			for (unsigned int i=2; i<4; i++)
			  isobar += evt->particle( i );
			//isobar += evt->particle(2);
			//isobar += evt->particle(4);

			/////////////////////////////////////////////////////////////////////
			// Do user calculation here for decay angles to match IsobarAngles //
			/////////////////////////////////////////////////////////////////////
			TLorentzVector beam = evt->particle(0);
			TLorentzVector recoil = evt->particle(1);
			TLorentzVector PX = resonance;
			TLorentzVector IsoBar2 = isobar;	
			TLorentzVector IsoBar1 = PX - IsoBar2;
			TLorentzVector IsoBar234 = IsoBar1 + evt->particle(3);
			TLorentzVector PBachX = PX - IsoBar2;
			
			// calculate decay angles in resonance X rest frame			
			TVector3 XRestBoost = PX.BoostVector();

			TLorentzVector beamX   = beam;
			TLorentzVector recoilX = recoil;
			TLorentzVector isobarX = IsoBar1;
			TLorentzVector isoX[4], bachX[4];;
			beamX.Boost(-1.0*XRestBoost);
			recoilX.Boost(-1.0*XRestBoost);
			isobarX.Boost(-1.0*XRestBoost);
			for(int i=0; i<4; i++) {
				bachX[i] = evt->particle(i+2);
				bachX[i].Boost(-1.0*XRestBoost);
				isoX[i] = PX - evt->particle(i+2);
				isoX[i].Boost(-1.0*XRestBoost);
			}

			// keep vectors for isobars in X rest frame for later angle definitions
			TLorentzVector temp, temp2, IsoBar1X, IsoBar2X, Pi0X, etaX, KpX, KmX;
			temp = IsoBar1; temp.Boost(-1.0*XRestBoost); IsoBar1X = temp;
			temp2 = IsoBar2; temp2.Boost(-1.0*XRestBoost); IsoBar2X = temp2;
			temp = evt->particle(2); temp.Boost(-1.0*XRestBoost); Pi0X = temp;
			temp = evt->particle(3); temp.Boost(-1.0*XRestBoost); etaX = temp;
			temp = evt->particle(4); temp.Boost(-1.0*XRestBoost); KpX = temp;
			temp = evt->particle(5); temp.Boost(-1.0*XRestBoost); KmX = temp;
			//---------------------------------------------------------------------

			// For GJ frame: choose beam as z-axis for reference 
			TVector3 z = beamX.Vect().Unit();
			TVector3 y = (beamX.Vect().Unit()).Cross((-recoilX.Vect().Unit())).Unit();
			TVector3 x = y.Cross(z);
	
			TVector3 anglesX( (isobarX.Vect()).Dot(x),
					  (isobarX.Vect()).Dot(y),
					  (isobarX.Vect()).Dot(z) );
		
			GDouble cosThetaBachX = anglesX.CosTheta();
			GDouble phiBachX = anglesX.Phi();

			GDouble cosThetaBachX_Iso[4];
			GDouble phiBachX_Iso[4];
			GDouble cosThetaBachIso1_Iso[4][3];
                        GDouble phiBachIso1_Iso[4][3];
			GDouble cosThetaBachIso2_Iso[4][3][2];
                        GDouble phiBachIso2_Iso[4][3][2];
			for(int i=0; i<4; i++) {
				TVector3 anglesX_Iso( (bachX[i].Vect()).Dot(x),
						      (bachX[i].Vect()).Dot(y),
						      (bachX[i].Vect()).Dot(z) );
				
				cosThetaBachX_Iso[i] = anglesX_Iso.CosTheta();
				phiBachX_Iso[i] = anglesX_Iso.Phi();

				// calculate decay angles in Isobar1 rest frame
				TVector3 Iso1RestBoost = isoX[i].BoostVector();	

				TLorentzVector isoIso1[3], bachIso1[3];
				int iso1count = 0;
                        	for(int j=0; j<4; j++) {
					if(i == j) continue; // skip bachelor from X decay
                                	bachIso1[iso1count] = evt->particle(j+2);
					bachIso1[iso1count].Boost(-1.0*XRestBoost);
					isoIso1[iso1count] = isoX[i] - bachIso1[iso1count];

					bachIso1[iso1count].Boost(-1.0*Iso1RestBoost);
                                        isoIso1[iso1count].Boost(-1.0*Iso1RestBoost);

					TVector3 zIso1 = isoX[i].Vect().Unit();
        		                TVector3 yIso1 = (z.Cross(zIso1)).Unit();
	                	        TVector3 xIso1 = yIso1.Cross(zIso1);

					TVector3 anglesIso1( (bachIso1[iso1count].Vect()).Dot(xIso1),
                                                             (bachIso1[iso1count].Vect()).Dot(yIso1),
                                                      	     (bachIso1[iso1count].Vect()).Dot(zIso1) );

	                                cosThetaBachIso1_Iso[i][iso1count] = anglesIso1.CosTheta();
        	                        phiBachIso1_Iso[i][iso1count] = anglesIso1.Phi();	
				

					// calculate decay angles in Isobar1 rest frame
	                                TVector3 Iso2RestBoost = isoIso1[iso1count].BoostVector();

        	                        TLorentzVector bachIso2[2];
					int iso2count = 0;
	                                for(int k=0; k<4; k++) {
                                        	if(i == k || j == k) continue; // skip bachelor from X and Iso1 decays
	                                        bachIso2[iso2count] = evt->particle(k+2);
        	                                bachIso2[iso2count].Boost(-1.0*XRestBoost);
                	                        bachIso2[iso2count].Boost(-1.0*Iso1RestBoost);
						bachIso2[iso2count].Boost(-1.0*Iso2RestBoost);
	
        	                                TVector3 zIso2 = isoIso1[iso1count].Vect().Unit();
                	                        TVector3 yIso2 = (zIso1.Cross(zIso2)).Unit();
                        	                TVector3 xIso2 = yIso2.Cross(zIso2);
	
        	                                TVector3 anglesIso2( (bachIso2[iso2count].Vect()).Dot(xIso2),
                	                                             (bachIso2[iso2count].Vect()).Dot(yIso2),
                        	                                     (bachIso2[iso2count].Vect()).Dot(zIso2) );
	
        	                                cosThetaBachIso2_Iso[i][iso1count][iso2count] = anglesIso2.CosTheta();
                	                        phiBachIso2_Iso[i][iso1count][iso2count] = anglesIso2.Phi();
						
						iso2count++;
					}

					iso1count++;
				}
			}

			//////////////////////////////////////
			// Test of successive 2-body decays //
			//////////////////////////////////////
			
/*
			TVector3 isoRestBoost1 = IsoBar1X.BoostVector();
			TLorentzVector PBachIso1 = KpX;
			TLorentzVector PBachIso1INV = KmX;
			PBachIso1.Boost(-1.0*isoRestBoost1);
			PBachIso1INV.Boost(-1.0*isoRestBoost1);
			
			TVector3 zIso1 = IsoBar1X.Vect().Unit();
			TVector3 yIso1 = (z.Cross(zIso1)).Unit();
			TVector3 xIso1 = yIso1.Cross(zIso1);

			TVector3 anglesIso1( (PBachIso1.Vect()).Dot(xIso1),
       					     (PBachIso1.Vect()).Dot(yIso1), 
					     (PBachIso1.Vect()).Dot(zIso1) );

			TVector3 anglesIso1INV( (PBachIso1INV.Vect()).Dot(xIso1),
						(PBachIso1INV.Vect()).Dot(yIso1), 
						(PBachIso1INV.Vect()).Dot(zIso1) );

			GDouble cosThetaBachIso1 = anglesIso1.CosTheta();
			GDouble phiBachIso1 = anglesIso1.Phi();
*/

			////////////////////////////////////////////////////////////
			// calculate decay angles in isobar rest (helicity) frame //
			////////////////////////////////////////////////////////////

			// KpKm restframe and decay
			TVector3 isoRestBoost1 = IsoBar1X.BoostVector();
			TLorentzVector PBachIso1 = KpX;
			TLorentzVector PBachIso1INV = KmX;
			PBachIso1.Boost(-1.0*isoRestBoost1);
			PBachIso1INV.Boost(-1.0*isoRestBoost1);
			
			TVector3 zIso1 = IsoBar1X.Vect().Unit();
			TVector3 yIso1 = (z.Cross(zIso1)).Unit();
			TVector3 xIso1 = yIso1.Cross(zIso1);

			TVector3 anglesIso1( (PBachIso1.Vect()).Dot(xIso1),
       					     (PBachIso1.Vect()).Dot(yIso1), 
					     (PBachIso1.Vect()).Dot(zIso1) );

			TVector3 anglesIso1INV( (PBachIso1INV.Vect()).Dot(xIso1),
						(PBachIso1INV.Vect()).Dot(yIso1), 
						(PBachIso1INV.Vect()).Dot(zIso1) );

			GDouble cosThetaBachIso1 = anglesIso1.CosTheta();
			GDouble phiBachIso1 = anglesIso1.Phi();
			
			GDouble cosThetaDiffIso1 = anglesIso1.Theta() - anglesIso1INV.Theta();
			GDouble phiDiffIso1 = anglesIso1.Phi() - anglesIso1INV.Phi();


			// Pi0eta restframe and decay
			TVector3 isoRestBoost2 = IsoBar2X.BoostVector();
			TLorentzVector PBachIso2 = Pi0X;
			TLorentzVector PBachIso2INV = etaX;
			PBachIso2.Boost(-1.0*isoRestBoost2);
			PBachIso2INV.Boost(-1.0*isoRestBoost2);
			
			TVector3 zIso2 = IsoBar2X.Vect().Unit();
			TVector3 yIso2 = (z.Cross(zIso2)).Unit();
			TVector3 xIso2 = yIso2.Cross(zIso2);

			TVector3 anglesIso2( (PBachIso2.Vect()).Dot(xIso2),
       					     (PBachIso2.Vect()).Dot(yIso2), 
					     (PBachIso2.Vect()).Dot(zIso2) );

			TVector3 anglesIso2INV( (PBachIso2INV.Vect()).Dot(xIso2),
						(PBachIso2INV.Vect()).Dot(yIso2), 
						(PBachIso2INV.Vect()).Dot(zIso2) );

			GDouble cosThetaBachIso2 = anglesIso2.CosTheta();
			GDouble phiBachIso2 = anglesIso2.Phi();
			
			GDouble cosThetaDiffIso2 = anglesIso2.Theta() - anglesIso2INV.Theta();
			GDouble phiDiffIso2 = anglesIso2.Phi() - anglesIso2INV.Phi();

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
					
					intenW->Fill( weightedInten );
					intenWVsM->Fill( resonance.M(), weightedInten );

					M_isobar->Fill( isobar.M() );
					M_recoil->Fill( recoil.M() );
					
					// calculate angular variables
					TLorentzVector beam = evt->particle ( 0 );
					TLorentzVector rec = evt->particle ( 1 );
					TLorentzVector p1 = evt->particle ( 2 );
					TLorentzVector target(0,0,0,rec[3]);
					
					if(isBaryonResonance) // assume t-channel
						t->Fill(-1*(beam-evt->particle(1)).M2());
					else
						t->Fill(-1*(recoil-target).M2());

					E->Fill(beam.E());
					EvsM->Fill(beam.E(),resonance.M());

					TLorentzRotation resonanceBoost( -resonance.BoostVector() );
					
					TLorentzVector beam_res = resonanceBoost * beam;
					TLorentzVector rec_res = resonanceBoost * rec;
					TLorentzVector p1_res = resonanceBoost * p1;
					
					// normal to the production plane
                                        TVector3 y = (beam.Vect().Unit().Cross(-rec.Vect().Unit())).Unit();

                                        // choose helicity frame: z-axis opposite recoil proton in rho rest frame
                                        TVector3 z = -1. * rec_res.Vect().Unit();
                                        TVector3 x = y.Cross(z).Unit();
                                        TVector3 angles( (p1_res.Vect()).Dot(x),
                                                         (p1_res.Vect()).Dot(y),
                                                         (p1_res.Vect()).Dot(z) );

                                        double cosTheta = angles.CosTheta();
                                        double phi = angles.Phi();

					M_CosTheta->Fill( resonance.M(), cosTheta);
					M_Phi->Fill( resonance.M(), phi);
					M_Phi_lab->Fill( resonance.M(), rec.Phi());
					
					TVector3 eps(1.0, 0.0, 0.0); // beam polarization vector
                                        double Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

					for(int i=0; i<4; i++) {
						M_CosTheta_Iso[i]->Fill( resonance.M(), cosThetaBachX_Iso[i]);
						M_Phi_Iso[i]->Fill( resonance.M(), phiBachX_Iso[i]);

						for(int j=0; j<3; j++) {
	                                                M_CosTheta_Iso1sq[i][j]->Fill( resonance.M(), cosThetaBachIso1_Iso[i][j]);
        	                                        M_Phi_Iso1sq[i][j]->Fill( resonance.M(), phiBachIso1_Iso[i][j]);

							for(int k=0; k<2; k++) {
                                                        	M_CosTheta_Iso2sq[i][j][k]->Fill( resonance.M(), cosThetaBachIso2_Iso[i][j][k]);
                                                        	M_Phi_Iso2sq[i][j][k]->Fill( resonance.M(), phiBachIso2_Iso[i][j][k]);
							}
						}
					}

					M_CosTheta_Iso1->Fill( IsoBar1.M(), cosThetaBachIso1);
					M_Phi_Iso1->Fill( IsoBar1.M(), phiBachIso1);
					M_Phi_lab_Iso1->Fill( IsoBar1.M(), IsoBar1.Phi());
					BeamE_CosTheta_Iso1->Fill(beam(3),cosThetaBachIso1);
					M_CosTheta_Diff_Iso1->Fill( IsoBar1.M(), cosThetaDiffIso1);
					M_Phi_Diff_Iso1->Fill( IsoBar1.M(), phiDiffIso1);
					ThetaCorr_Iso1->Fill(anglesIso1.Theta(), anglesIso1INV.Theta());
					PhiCorr_Iso1->Fill(anglesIso1.Phi(), anglesIso1INV.Phi());

					M_CosTheta_Iso2->Fill( IsoBar2.M(), cosThetaBachIso2);
					M_Phi_Iso2->Fill( IsoBar2.M(), phiBachIso2);
					M_Phi_lab_Iso2->Fill( IsoBar2.M(), IsoBar2.Phi());
					M_CosTheta_Diff_Iso2->Fill( IsoBar2.M(), cosThetaDiffIso2);
					M_Phi_Diff_Iso2->Fill( IsoBar2.M(), phiDiffIso2);
					ThetaCorr_Iso2->Fill(anglesIso2.Theta(), anglesIso2INV.Theta());
					PhiCorr_Iso2->Fill(anglesIso2.Phi(), anglesIso2INV.Phi());

					for(int ipart=1; ipart<6; ipart++) 
						PVsTheta[ipart]->Fill(evt->particle(ipart).Theta(), evt->particle(ipart).Vect().Mag());
					
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
				TLorentzVector rec = evt->particle ( 1 );
				
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

	for(int i=0; i<4; i++) {
		M_CosTheta_Iso[i]->Write();
		M_Phi_Iso[i]->Write();
		
		for(int j=0; j<3; j++) {
         	       M_CosTheta_Iso1sq[i][j]->Write();
                       M_Phi_Iso1sq[i][j]->Write();
			
		       for(int k=0; k<2; k++) {
                      	     M_CosTheta_Iso2sq[i][j][k]->Write();
                       	     M_Phi_Iso2sq[i][j][k]->Write();
		       }
	        }
	}

	M_CosTheta_Iso1->Write();
	M_Phi_Iso1->Write();
	M_Phi_lab_Iso1->Write();
	M_CosTheta_Diff_Iso1->Write();
	M_Phi_Diff_Iso1->Write();
	M_CosTheta_Iso2->Write();
	M_Phi_Iso2->Write();
	M_Phi_lab_Iso2->Write();
	BeamE_CosTheta_Iso1->Write();
	M_CosTheta_Diff_Iso2->Write();
	M_Phi_Diff_Iso2->Write();

	ThetaCorr_Iso1->Write();
	PhiCorr_Iso1->Write();
	ThetaCorr_Iso2->Write();
	PhiCorr_Iso2->Write();

	for(int ipart=1; ipart<6; ipart++) 
		PVsTheta[ipart]->Write();

	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}

