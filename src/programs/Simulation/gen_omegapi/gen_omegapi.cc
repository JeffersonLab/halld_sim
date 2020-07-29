//Generator for b1(1235)->omega pi by A. M. Foda https://github.com/amfodajlab
//Based on gen_amp by Alex Austregesilo https://github.com/aaust
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
#include "AMPTOOLS_DATAIO/ASCIIDataWriter.h"

#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/omegapiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToNPartP.h"
#include "AMPTOOLS_MCGEN/NBodyPhaseSpaceFactory.h"

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

#include "UTILITIES/CobremsGeneration.hh"
#include "UTILITIES/BeamProperties.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
  
	string  configfile("");
	string  outname("");
	string  hddmname("");
        string  asciiname("");
	
	bool diag = false;
	bool genFlat = false;
	
	// default upper and lower bounds 
	double lowMass = 1.0;//To take over threshold with a BW omega mass
	double highMass = 2.0;

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
	
	float Mpip=ParticleMass(PiPlus), Mpi0=ParticleMass(Pi0), Momega=0.782;

	//Exprected particle list: 
	// pi0 omega(pi0 "rho"(pi+ pi-))
	//  2         3         4   5
	int par_types_list[]={1,14,7,7,8,9};
	vector<int> part_types(par_types_list,par_types_list+6);

	float part_masses_list1[]={Mpi0, Momega};
	vector<double> part_masses1(part_masses_list1,part_masses_list1+2);

	float part_masses_list2[]={Mpip, Mpip, Mpi0};
	vector<double> part_masses2(part_masses_list2,part_masses_list2+3);

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
                if (arg == "-oascii"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  asciiname = argv[++i]; }
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
		cout << "No config file or output specificed:  run gen_omegapi -h for help" << endl;
		exit(1);
	}
	
	// open config file and be sure only one reaction is specified
	ConfigFileParser parser( configfile );
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	assert( cfgInfo->reactionList().size() == 1 );
	ReactionInfo* reaction = cfgInfo->reactionList()[0];

	// check for unstable particle at lower vertex
        vector<Particle_t> ParticlesLowerVertex;
        vector<double> massesLowerVertex;
        double thresholdLowerVertex = 0;
        vector< BreitWignerGenerator > bwGenLowerVertex;
        const vector<ConfigFileLine> configFileLinesLowerVertex = parser.getConfigFileLines();
        for (vector<ConfigFileLine>::const_iterator it=configFileLinesLowerVertex.begin(); it!=configFileLinesLowerVertex.end(); it++) {
          if ((*it).keyword() == "define"){
            if( (*it).arguments()[0] != "lowerVertex") continue;
            bwGenLowerVertex.push_back( BreitWignerGenerator( atof((*it).arguments()[1].c_str()), atof((*it).arguments()[2].c_str())) );
            cout << "Unstable particle at lower vertex: mass = " << (*it).arguments()[1]  << "GeV , width = " << (*it).arguments()[2] << "GeV" << endl;
            for(unsigned int i=3; i<(*it).arguments().size(); i++) {
              ParticlesLowerVertex.push_back(ParticleEnum((*it).arguments()[i].c_str()));
              massesLowerVertex.push_back(ParticleMass(ParticlesLowerVertex[i-3]));
              thresholdLowerVertex += ParticleMass(ParticlesLowerVertex[i-3]);
            }
            break;
	  }
	}
	
	// use particletype.h to convert reaction particle names (for upper vertex)
	vector<Particle_t> Particles;
	// don't include non-nucleon lower vertex decay particles in meson decay
	unsigned int maxUpperVertexChild = reaction->particleList().size();
        if(bwGenLowerVertex.size() == 1) maxUpperVertexChild -= (ParticlesLowerVertex.size()-1);
        for (unsigned int i = 0; i < maxUpperVertexChild; i++){
	  Particle_t locEnum = ParticleEnum(reaction->particleList()[i].c_str());
	  // Beam particle is always photon
	  if (locEnum == 0 && i > 0)
	    cout << "ConfigFileParser WARNING:  unknown particle type \"" << reaction->particleList()[i] << "\"" << endl;
	  Particles.push_back(ParticleEnum(reaction->particleList()[i].c_str()));
      }

	vector<double> childMasses;
	double threshold = 0;
	if(bwGenLowerVertex.size() == 0) childMasses.push_back(0.135); // b1 -> omega pi0
	else childMasses.push_back(0.1396); // b1 -> omega pi+/-
	childMasses.push_back(0.782);
	threshold = 0.135 + 0.782;

	// loop to look for resonance in config file
	// currently only one at a time is supported 
	const vector<ConfigFileLine> configFileLines = parser.getConfigFileLines();
	double resonance[]={1.235, 0.142};
	bool foundResonance = false;

	for (vector<ConfigFileLine>::const_iterator it=configFileLines.begin(); it!=configFileLines.end(); it++) {
	  if ((*it).keyword() == "define") {
	    if ((*it).arguments()[0] == "rho" || (*it).arguments()[0] == "omega" || (*it).arguments()[0] == "phi" || (*it).arguments()[0] == "b1" || (*it).arguments()[0] == "a1" || (*it).arguments()[0] == "Lambda1520"){
	      if ( (*it).arguments().size() != 3 )
		continue;
	      resonance[0]=atof((*it).arguments()[1].c_str());
	      resonance[1]=atof((*it).arguments()[2].c_str());
	      cout << "Distribution seeded with resonance " << (*it).arguments()[0] << " : mass = " << resonance[0] << "GeV , width = " << resonance[1] << "GeV" << endl; 
	      if((*it).arguments()[0] == "Lambda1520")
	      foundResonance = true;
	      break;
	    }
	  }
	}
	if (!foundResonance)
	  cout << "ConfigFileParser WARNING:  no known resonance found, seed with mass = 1.235, width = 0.142 GeV" << endl; 

	// random number initialization (set to 0 by default)
	TRandom3* gRandom = new TRandom3();
	gRandom->SetSeed(seed);
	seed = gRandom->GetSeed();
	cout << "TRandom3 Seed : " << seed << endl;

	// setup AmpToolsInterface
	AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
        AmpToolsInterface::registerAmplitude( BreitWigner() );
        AmpToolsInterface::registerAmplitude( Uniform() );


	AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );

	double polAngle = -1;//amorphous
	// loop to look for beam configuration file
        TString beamConfigFile;
        const vector<ConfigFileLine> configFileLinesBeam = parser.getConfigFileLines();
        for (vector<ConfigFileLine>::const_iterator it=configFileLinesBeam.begin(); it!=configFileLinesBeam.end(); it++) {
                if ((*it).keyword() == "define") {
                        TString beamArgument =  (*it).arguments()[0].c_str();
                        if(beamArgument.Contains("beamconfig")) {
                                beamConfigFile = (*it).arguments()[1].c_str();
		  		BeamProperties beamProp(beamConfigFile);
				polAngle = beamProp.GetPolAngle();
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
	GammaPToNPartP resProd = GammaPToNPartP( threshold<lowMass ? lowMass : threshold, highMass, childMasses, ProductionMechanism::kProton, type, slope, lowT, highT, seed, beamConfigFile );

	vector< BreitWignerGenerator > m_bwGen;
        m_bwGen.push_back( BreitWignerGenerator(0.782, 0.008) );
		
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
	
	ASCIIDataWriter* asciiOut = NULL;
        if( asciiname.size() != 0 ) asciiOut = new ASCIIDataWriter( asciiname );

	TFile* diagOut = new TFile( "gen_omegapi_diagnostic.root", "recreate" );
	ostringstream locStream;
	ostringstream locIsobarStream;
	ostringstream locIsobar2Stream;
	for (unsigned int i=2; i<Particles.size(); i++){
	  locStream << ParticleName_ROOT(Particles[i]);
	  if ( i> 2 )
	    locIsobarStream << ParticleName_ROOT(Particles[i]);
	}
	string locHistTitle = string("Resonance Mass ;") + locStream.str() + string(" Invariant Mass (GeV/c^{2});");
	string locIsobarTitle = string("Isobar Mass ;") + locIsobarStream.str() + string(" Invariant Mass (GeV/c^{2});");
	string locIsobar2Title = string("Isobar2 Mass ;") + locIsobar2Stream.str() + string(" Invariant Mass (GeV/c^{2});");

	TH1F* mass = new TH1F( "M", locHistTitle.c_str(), 180, lowMass, highMass );
	TH1F* massW = new TH1F( "M_W", ("Weighted "+locHistTitle).c_str(), 180, lowMass, highMass );
	massW->Sumw2();
	TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
	TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
	
	TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );

	TH1F* M_isobar = new TH1F( "M_isobar", locIsobarTitle.c_str(), 200, 0, 2 );
	TH1F* M_isobar2 = new TH1F( "M_isobar2", locIsobar2Title.c_str(), 200, 0, 2 );
	TH1F* M_recoil = new TH1F( "M_recoil", "; Recoil mass (GeV)", 200, 0, 2 );
	TH1F* M_p1 = new TH1F( "M_p1", "p1", 200, 0, 2 );
	TH1F* M_p2 = new TH1F( "M_p2", "p2", 200, 0, 2 );
	TH1F* M_p3 = new TH1F( "M_p3", "p3", 200, 0, 2 );
	TH1F* M_p4 = new TH1F( "M_p4", "p4", 200, 0, 2 );

        TH2F* M_dalitz = new TH2F( "M_dalitz", "dalitzxy", 200, -5, 5, 200, -5, 5);

	TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
	TH2F* M_CosTheta = new TH2F( "M_CosTheta", "M vs. cos#vartheta", 180, lowMass, highMass, 200, -1, 1);
	TH2F* M_Phi = new TH2F( "M_Phi", "M vs. #varphi", 180, lowMass, highMass, 200, -3.14, 3.14);
	TH2F* M_CosThetaH = new TH2F( "M_CosThetaH", "M vs. cos#vartheta_{H}", 180, lowMass, highMass, 200, -1, 1);
	TH2F* M_PhiH = new TH2F( "M_PhiH", "M vs. #varphi_{H}", 180, lowMass, highMass, 200, -3.14, 3.14);
	TH2F* M_Phi_Prod = new TH2F( "M_Phi_Prod", "M vs. #Phi_{Prod}", 180, lowMass, highMass, 200, -3.14, 3.14);

	int eventCounter = 0;
	while( eventCounter < nEvents ){
		
		if( batchSize < 1E4 ){
			
			cout << "WARNING:  small batches could have batch-to-batch variations\n"
			     << "          due to different maximum intensities!" << endl;
		}
		
		cout << "Generating four-vectors..." << endl;

		// decay omega (and Delta++, if generated)
		ati.clearEvents();
		for( int i = 0; i < batchSize; ++i ){

			// setup omega decay
                        pair< double, double > bw = m_bwGen[0]();
                        double omega_mass_bw = bw.first;
                        if ( omega_mass_bw < 0.45 || omega_mass_bw > 0.864) continue;//Avoids Tcm < 0 in NBPhaseSpaceFactory and BWgenerator

			vector<double> childMasses_omega_bw;
              		childMasses_omega_bw.push_back(childMasses[0]);
        		childMasses_omega_bw.push_back(omega_mass_bw);
			
			resProd.setChildMasses(childMasses_omega_bw);
			resProd.getProductionMechanism().setMassRange( lowMass, highMass );

			// setup lower vertex decay
			pair< double, double > bwLowerVertex;
			double lowerVertex_mass_bw = 0.;
			if(bwGenLowerVertex.size() == 1) {
				bwLowerVertex = bwGenLowerVertex[0]();
				lowerVertex_mass_bw = bwLowerVertex.first;
				if ( lowerVertex_mass_bw < thresholdLowerVertex || lowerVertex_mass_bw > 2.0) continue;
				resProd.getProductionMechanism().setRecoilMass( lowerVertex_mass_bw );
			}

			  Kinematics* step1 = resProd.generate();
			  TLorentzVector beam = step1->particle( 0 );
			  TLorentzVector recoil = step1->particle( 1 );
			  TLorentzVector bachelor_pi = step1->particle( 2 );
			  TLorentzVector omega = step1->particle( 3 );
			  TLorentzVector b1 = bachelor_pi + omega;

			  // decay step for omega
			  //cout << "omega mass =" << omega_mass_bw << endl;
        		  NBodyPhaseSpaceFactory omega_to_pions = NBodyPhaseSpaceFactory( omega_mass_bw, part_masses2);
			  vector<TLorentzVector> omega_daughters = omega_to_pions.generateDecay();
			  
			  TLorentzVector piplus = omega_daughters[0];//second decay step
			  TLorentzVector piminus = omega_daughters[1];//second decay step
			  TLorentzVector omegas_pi0 = omega_daughters[2];//second decay step
			  //cout << "second step masses ="<< piplus.M() << ", "<< piminus.M() << ", " << omegas_pi0.M() << endl;

			  omegas_pi0.Boost( omega.BoostVector() );
			  piplus.Boost( omega.BoostVector() );			  
			  piminus.Boost( omega.BoostVector() );
		  
			  // decay step for Delta++
			  TLorentzVector nucleon;
			  vector<TLorentzVector> lowerVertexChild;
			  if(bwGenLowerVertex.size() == 1) {
				  NBodyPhaseSpaceFactory lowerVertex_decay = NBodyPhaseSpaceFactory( lowerVertex_mass_bw, massesLowerVertex);
				  lowerVertexChild = lowerVertex_decay.generateDecay();

				  // boost to lab frame via recoil kinematics
				  for(unsigned int j=0; j<lowerVertexChild.size(); j++)
					  lowerVertexChild[j].Boost( recoil.BoostVector() );
				  nucleon = lowerVertexChild[0];	
			  }		  
			  else 
				  nucleon = recoil;

			  // store particles in kinematic class
			  vector< TLorentzVector > allPart;
			  //same order as config file, omegapi Amplitudes and ReactionFilter
			  allPart.push_back( beam );
			  allPart.push_back( nucleon );
			  allPart.push_back( bachelor_pi );
			  allPart.push_back( omegas_pi0 );
			  allPart.push_back( piplus );
			  allPart.push_back( piminus );
			  if(bwGenLowerVertex.size() == 1)
				  for(unsigned int j=1; j<lowerVertexChild.size(); j++)
                                        allPart.push_back(lowerVertexChild[j]);

			  Kinematics* kin = new Kinematics( allPart, 1.0 );
			  ati.loadEvent( kin, i, batchSize );
			  delete step1;
			  delete kin;
    		}
		
		cout << "Processing events..." << endl;
		
		// include factor of 1.5 to be safe in case we miss peak -- avoid
		// intensity calculation of we are generating flat data
		double maxInten = ( genFlat ? 1 : 1.50* ati.processEvents( reaction->reactionName() ) );
		
		
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* evt = ati.kinematics( i );
			TLorentzVector resonance;
			for (unsigned int i=2; i<Particles.size(); i++)
			  resonance += evt->particle( i );

			TLorentzVector isobar;
			for (unsigned int i=3; i<Particles.size(); i++)
			  isobar += evt->particle( i );
			
			TLorentzVector isobar2;
			for (unsigned int i=4; i<Particles.size(); i++)
			  isobar2 += evt->particle( i );
			
			TLorentzVector recoil = evt->particle( 1 );
                        if(bwGenLowerVertex.size()) {
				for(unsigned int j=Particles.size(); j<evt->particleList().size(); j++)
                                        recoil += evt->particle( j );
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
					
					intenW->Fill( weightedInten );
					intenWVsM->Fill( resonance.M(), weightedInten );

					M_isobar->Fill( isobar.M() );
					M_isobar2->Fill( isobar2.M() );
					M_recoil->Fill( recoil.M() );
					
					// calculate angular variables
					TLorentzVector beam = evt->particle ( 0 );
					TLorentzVector p1 = evt->particle ( 2 );
					TLorentzVector p2 = evt->particle ( 3 );
					TLorentzVector p3 = evt->particle ( 4 );
					TLorentzVector p4 = evt->particle ( 5 );
					TLorentzVector target(0,0,0,ParticleMass(Proton));
					
					M_p1->Fill( p1.M() );
					M_p2->Fill( p2.M() );
					M_p3->Fill( p3.M() );
					M_p4->Fill( p4.M() );

					double dalitz_s, dalitz_t, dalitz_u, dalitz_d, dalitz_sc, dalitzx, dalitzy;
					dalitz_s = (p3+p2).M2();//s=M(pip pi0)
					dalitz_t = (p4+p2).M2();//s=M(pim pi0)
					dalitz_u = (p3+p4).M2();//s=M(pip pim)
					dalitz_d = 2*(p2+p3+p4).M()*( (p2+p3+p4).M() - ((2*139.57018)+134.9766) );
					dalitz_sc = (1/3)*( (p2+p3+p4).M2() - ((2*(139.57018*139.57018))+(134.9766*134.9766)) );
					dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
					dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;

					M_dalitz->Fill(dalitzx,dalitzy);
					
					t->Fill(-1*(recoil-target).M2());

                                        TLorentzVector Gammap = beam + target;
                                        vector <double> loccosthetaphi = getomegapiAngles(polAngle, isobar, resonance, beam, Gammap);
                                        double cosTheta = cos(loccosthetaphi[0]);
                                        double phi = loccosthetaphi[1];

                                        vector <double> loccosthetaphih = getomegapiAngles( p3, isobar, resonance, Gammap, p4);
                                        double cosThetaH = cos(loccosthetaphih[0]);
                                        double phiH = loccosthetaphih[1];

					M_CosTheta->Fill( resonance.M(), cosTheta);
					M_Phi->Fill( resonance.M(), phi);
					M_CosThetaH->Fill( resonance.M(), cosThetaH);
					M_PhiH->Fill( resonance.M(), phiH);

                                        double Phi = loccosthetaphi[2];
					M_Phi_Prod->Fill( resonance.M(), Phi);

                                        GDouble psi = phi - Phi;
                                        if(psi < -1*PI) psi += 2*PI;
                                        if(psi > PI) psi -= 2*PI;
					
					CosTheta_psi->Fill( psi, cosTheta);
										
					// we want to save events with weight 1
					evt->setWeight( 1.0 );
					
					if( hddmOut ) hddmOut->writeEvent( *evt, pTypes );
                                        if( asciiOut ) asciiOut->writeEvent( *evt, pTypes );
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
	M_isobar2->Write();
	M_recoil->Write();
        M_p1->Write();
        M_p2->Write();
        M_p3->Write();
        M_p4->Write();
        M_dalitz->Write();
	t->Write();
	CosTheta_psi->Write();
	M_CosTheta->Write();
	M_Phi->Write();
	M_CosThetaH->Write();
	M_PhiH->Write();
	M_Phi_Prod->Write();

	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}
