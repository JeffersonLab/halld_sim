
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

#include "AMPTOOLS_DATAIO/DataWriter.h"
#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "AMPTOOLS_DATAIO/FSRootDataWriter.h"
#include "AMPTOOLS_MCGEN/HDDMDataWriter.h"

#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/BreitWigner3body.h"
#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
#include "AMPTOOLS_AMPS/VecRadiative_SDME.h"
#include "AMPTOOLS_AMPS/Lambda1520Angles.h"
#include "AMPTOOLS_AMPS/Lambda1520tdist.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/Ylm.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/DblRegge_FastEta.h"
#include "AMPTOOLS_AMPS/DblRegge_FastPi.h"
#include "AMPTOOLS_AMPS/Hist2D.h"
#include "AMPTOOLS_AMPS/Flatte.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"

#include "AMPTOOLS_MCGEN/FixedTargetGenerator.h"
#include "AMPTOOLS_MCGEN/BreitWignerGenerator.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

#include "UTILITIES/BeamProperties.h"


#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

using std::complex;
using namespace std;


vector<int> parseString(string vertexString);
string checkParticle(string particleString);
inline static unsigned short int UseResonance(Particle_t p);
int main( int argc, char* argv[] ){
  
  TString beamConfigFile("");
  string  atConfigFile("");
  string  hddmName("");
  string  outName("");
  string  lvString("");
  string  uvString("");

  	
  bool fixedGen = false;
  bool fsRootFormat = false;
  bool diag = false;

  // default upper and lower bounds -- these
  // just need to be outside the kinematic bounds
  // of any reaction
  
  double beamMaxE   = 12.0;
  double beamPeakE  = 9.0;
  double beamLowE   = 3.0;
  double beamHighE  = 12.0;

  int runNum = 30731;
  unsigned int seed = 0;
  unsigned int reWeight =7;

  int nEvents = 10000;
  int batchSize = 10000;
  
  vector<int> indicateBW;
  map<int,BreitWignerGenerator> mpBW;
  mpBW[PDGtype( ParticleEnum("Omega") )] = BreitWignerGenerator( ParticleMass( ParticleEnum( "Omega" ) ), 0.00868, seed); // Initialize BW for omega
  mpBW[PDGtype( ParticleEnum("Phi") )] = BreitWignerGenerator( ParticleMass( ParticleEnum( "Phi" ) ), 0.00868, seed); // Initialize BW for omega

  // Initialization of FixedTargetGenerator
  FixedTargetGenerator ftGen;

  //parse command line:
  for (int i = 1; i < argc; i++){
		
    string arg(argv[i]);
		
    if (arg == "-ac"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  atConfigFile = argv[++i]; }
    if (arg == "-bc"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  beamConfigFile = argv[++i]; }
    if (arg == "-o"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outName = argv[++i]; }
    if (arg == "-uvBW") {
      if ((i+1 == argc) || (i+2 == argc) || (i+3 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-') || (argv[i+3][0] == '-')) arg = "-h";
      else{ 
	ftGen.addUpperVtxBW( atof( argv[i+1] ), atof( argv[i+2] ), atof( argv[i+3] ) ); 
      }
    }
    if (arg == "-lvBW") {
      if ((i+1 == argc) || (i+2 == argc) || (i+3 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-') || (argv[i+3][0] == '-')) arg = "-h";
      else{ 
        ftGen.addLowerVtxBW( atof( argv[i+1] ), atof( argv[i+2] ), atof( argv[i+3] ) );	
      }
    }
    if (arg == "-uv"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  uvString = argv[++i]; }
    if (arg == "-lv"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  lvString = argv[++i]; }
    if (arg == "-fsroot") fsRootFormat = true; 
    if (arg == "-hd"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  hddmName = argv[++i]; }
    if (arg == "-n"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  nEvents = atoi( argv[++i] ); }
    if (arg == "-r"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  runNum = atoi( argv[++i] ); }
    if (arg == "-s"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  seed = atoi( argv[++i] ); }
    if (arg == "-lvRange"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
	ftGen.setLowerVtxRange( atof( argv[i+1] ), atof( argv[i+2] ) ); 
      }
    }
    if (arg == "-uvRange"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
	ftGen.setUpperVtxRange( atof( argv[i+1] ), atof( argv[i+2] ) ); 
      }
    }
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
    if (arg == "-t"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
        ftGen.addMomentumTransfer( atof( argv[i+1] ), atof( argv[i+2] ) );
      }
    }
    if (arg == "-tRange"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
        ftGen.setMomentumTransferRange( atof( argv[i+1] ), atof( argv[i+2] ) );
      }
    }
    if (arg == "-mask"){
      if ((i+1 == argc) || (i+2 == argc) || (i+3 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-') || (argv[i+3][0] == '-')) arg = "-h";
      else{
      	if( atoi( argv[i+1] ) == 0 ) reWeight -= FixedTargetGenerator::kUpperVtxMass;	
      	if( atoi( argv[i+2] ) == 0 ) reWeight -= FixedTargetGenerator::kLowerVtxMass;
	if( atoi( argv[i+3] ) == 0 ) reWeight -= FixedTargetGenerator::kMomentumTransfer;
      }
    }
    if (arg == "-f") fixedGen = true; 
    if (arg == "-d") diag = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -ac       <file>" << endl;
      cout << "\t           AmpTools config file" << endl;
      cout << "\t -bc       <file>" << endl;
      cout << "\t           beam config file (will generate local file within coherent peak if not included) " << endl;
      cout << "\t -o        <name>" << endl;
      cout << "\t           ROOT file output name" << endl;
      cout << "\t -uv       <string>" << endl;
      cout << "\t           upper-vertex indices (0 is first) in AmpTools reaction" << endl;
      cout << "\t -lv       <string>" << endl;
      cout << "\t           lower-vertex indices (0 is first) in AmpTools reaction" << endl;
      cout << "\t -fsroot"         << endl;
      cout << "\t           Enable output in FSRoot format [optional]" << endl;
      cout << "\t -hd       <name> " << endl;
      cout << "\t           HDDM file output name [optional]" << endl;
      cout << "\t -f " << endl;
      cout << "\t           Option will skip AmpTools accept/reject process. Output information of FixedTargetGenerator [optional]" << endl;
      cout << "\t -d "            << endl; 
      cout << "\t           Plot only diagnostic histograms [optional]" << endl; 
      cout << "\t -n        <number>" << endl;
      cout << "\t           Minimum number of events to generate [optional]" << endl;
      cout << "\t -r        <run number>" << endl;
      cout << "\t           Run number assigned to generated events [optional]" << endl;
      cout << "\t -s        <number>" << endl;
      cout << "\t           Random number seed initialization [optional]" << endl;
      cout << "\t -lvRange  <min> <max> " << endl;
      cout << "\t           Lower vertex mass range ( lowerLimit upperLimit ) (GeV) [optional]" << endl;
      cout << "\t -uvRange  <min> <max> " << endl; 
      cout << "\t           Upper vertex mass range  ( lowerLimit upperLimit ) (GeV) [optional]" << endl;
      cout << "\t -t        <tSlope> <scale>" << endl;
      cout << "\t           Momentum transfer slope ( t Scale ) [optional]" << endl;
      cout << "\t -tRange   <min> <max>" << endl;
      cout << "\t           Momentum transfer range (lowerLimit upperLimit) (Gev)^2 [optional]" << endl;
      cout << "\t -uvBW     <mass> <width> <scale> " << endl;
      cout << "\t           Concentrate events around one or more upper vertex resonances (Mass Width Scale) [optional]" << endl;
      cout << "\t -lvBW     <mass> <width> <scale> " << endl;
      cout << "\t           Concentrate events around one or more lower vertex resonances (Mass Width Scale) [optional]" << endl;
      cout << "\t -mask     <upper vertex mass> <lower vertex mass> <t distribution>" << endl;
      cout << "\t           Remove weight to recover phase space for upper vertex mass, lower vertex, or momentum transfer t respectively ( uvM lvM t ) [optional]" << endl; 
      cout << "\t -m        <energy> " << endl;
      cout << "\t           Electron beam energy (or photon energy endpoint) [optional]" << endl;
      cout << "\t -p        <peak> " << endl;
      cout << "\t           Coherent peak photon energy [optional]" << endl;
      cout << "\t -a        <min E beam>" << endl;
      cout << "\t           Minimum photon energy to simulate events [optional]" << endl;
      cout << "\t -b        <max E beam>" << endl;
      cout << "\t           Maximum photon energy to simulate events [optional]" << endl << endl;
      exit(1);
    }
  }
 	
  if( atConfigFile.size() == 0 || outName.size() == 0 ||
      lvString.size() == 0 || uvString.size() == 0 ){
    
    cout << "ERROR:  Missing required arguments:  run gen_amp_v2 -h for help" << endl;
    exit( 1 );
  }
	
  // open config file and be sure only one reaction is specified
  ConfigFileParser parser( atConfigFile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();

  // the logic that follows assumes there is only one reaction defined
  assert( cfgInfo->reactionList().size() == 1 );
  ReactionInfo* reaction = cfgInfo->reactionList()[0];

  //

  vector< int > lvIndices = parseString( lvString );
  vector< int > uvIndices = parseString( uvString );
  vector< int > valPDG;
  vector< int > valVertex;
  // get the masses for each.....

  vector< double > lvMasses;
  vector< double > uvMasses;
  vector< int > pTypes;
  ostringstream upperVtxName;
  ostringstream lowerVtxName;
 
  for( unsigned int i = 0; i < reaction->particleList().size(); ++i){
    // checkParticle will check if an identical particle is used, dictated by % after particle name
    // if none then will return string unchanged
    string tempString = checkParticle( reaction->particleList()[i] );
    Particle_t particle = ParticleEnum( tempString.c_str() );
    if (particle == 0 && i > 0){
      cout << "ERROR:  unknown particle " << tempString 
           << " unable to configure generator." << endl;
      exit( 1 );
    }
    else{
      // this loop will check if particle is part of lower vertex according to user
      for( unsigned int j = 0; j < lvIndices.size(); ++j){
        if( lvIndices[j] > 0 && ((unsigned int)lvIndices[j] == i )){
       	  cout << "This is lower vertex particle with indices " << lvIndices[j] 
	  << " with name " <<  tempString << endl; 
	  lvMasses.push_back( ParticleMass( particle ) );	
	  lowerVtxName << ParticleName_ROOT( particle );
	  if( UseResonance( particle ) ) {   // Check if particle has a resonance 
            cout << "Particle with resonance found in lower vertex!" << endl << endl;
	    valPDG.push_back( PDGtype( particle ) );
            valVertex.push_back( (int)lvMasses.size() - 1 );
            indicateBW.push_back(1); //This will vary lower BW mass
          }
        }
      }
      // this loop will check if particle is part of upper vertex according to user
      for( unsigned int k = 0; k < uvIndices.size(); ++k){
        if( uvIndices[k] > 0 && ((unsigned int)uvIndices[k] == i )){
          cout << "This is upper vertex particle with indices " << uvIndices[k] 
          << " with name " <<  tempString << endl;
	  uvMasses.push_back( ParticleMass( particle ) ); 
          upperVtxName << ParticleName_ROOT( particle );  
	  if( UseResonance( particle ) ){
            cout << "Particle with resonance found in upper vertex!" << endl << endl;
            valPDG.push_back( PDGtype( particle ) );
            valVertex.push_back( (int)uvMasses.size() - 1 );
            indicateBW.push_back(2); //This will vary upper BW mass
          }   
        }
      }
      // add particle to pTypes
      pTypes.push_back( particle );
    }
    
  }
  map<int, BreitWignerGenerator>::iterator itBW = mpBW.begin();
 

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
  AmpToolsInterface::registerAmplitude( VecRadiative_SDME() );
  AmpToolsInterface::registerAmplitude( Lambda1520Angles() );
  AmpToolsInterface::registerAmplitude( Lambda1520tdist() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
  AmpToolsInterface::registerAmplitude( Ylm() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( Hist2D() );
  AmpToolsInterface::registerAmplitude( DblRegge_FastEta() );
  AmpToolsInterface::registerAmplitude( DblRegge_FastPi() );
  AmpToolsInterface::registerAmplitude( Flatte() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( ComplexCoeff() );
  AmpToolsInterface::registerAmplitude( OmegaDalitz() );
  AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );

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
	
  // get beam properties from configuration file
  BeamProperties beamProp( beamConfigFile );
  TH1D* cobrem_vs_E = (TH1D*)beamProp.GetFlux();

  // this will need reset every event, but start with a value
  // to setup the generator
  double beamEnergy = cobrem_vs_E->GetRandom();

  double targetMass = ParticleMass( ParticleEnum( "Proton" ) );
  ftGen.setTargetMass( targetMass );
  ftGen.setBeamEnergy( beamEnergy );
  ftGen.setUpperVtxMasses( uvMasses );
  ftGen.setLowerVtxMasses( lvMasses );
  ftGen.setSeed( seed );
  // Add the new seed value to sub BW's if exist
  mpBW[PDGtype( ParticleEnum("omega") )].setSeed( seed );
  mpBW[PDGtype( ParticleEnum("phi") )].setSeed( seed );
 
  // Sets reweighting based off of options given in command line
  ftGen.setReweightMask( reWeight );
  
  HDDMDataWriter* hddmOut = NULL;
  if( hddmName.size() != 0 ) hddmOut = new HDDMDataWriter( hddmName, runNum, seed);
  
  // the first argument to the FSRootDataWriter is the number of particles *in addition to* the beam
  // particle, which is typically the first in the list in GlueX reaction definitions
  DataWriter* rootOut = ( fsRootFormat ?
			  static_cast< DataWriter*>( new FSRootDataWriter( reaction->particleList().size()-1, outName ) ) :
			  static_cast< DataWriter* >( new ROOTDataWriter( outName ) ) );
	
  TFile* diagOut = new TFile( "gen_amp_diagnostic.root", "recreate" );
  
  string locHistTitle = string("Meson Mass ;") + upperVtxName.str() + string(" Invariant Mass (GeV/c^{2});");
  string locRecoilTitle = string("Baryon Mass ;") + lowerVtxName.str() + string(" Invariant Mass (GeV/c^{2});");

	
  TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );
  TH1F* tW = new TH1F( "tW", "-t Distribution", 200, 0, 2 ); 
  TH1F* e = new TH1F( "e", "Beam Energy", 120, 0, 12 );
  TH1F* eW = new TH1F( "eW", "Beam Energy", 120, 0, 12 );
  TH1F* eWI = new TH1F( "eWI", "Beam Energy", 120, 0, 12 );
  TH1F* intenW = new TH1F("intenW", "True PDF/ Gen. PDF", 1000, -0.1e-03, 0.8e-03);
  TH2F* intenWVsE = new TH2F("intenWVsE","Ratio vs. E", 100, 0, 12, 200, -0.1e-03, 0.8e-03);
  TH2F* intenWVsM = new TH2F("intenWVsE","Ratio vs. M", 100, 0, 3., 200, -0.1e-03, 0.8e-03);
  
  TH1F* m_Meson = new TH1F( "m_Meson", locHistTitle.c_str(), 200, 0., 3. );
  TH1F* mW_Meson = new TH1F( "mW_Meson", locHistTitle.c_str(), 200, 0., 3. );
  TH1F* m_recoil = new TH1F( "m_recoil", locRecoilTitle.c_str(), 200, 0., 3. );
  TH1F* mW_recoil = new TH1F( "mW_recoil", locRecoilTitle.c_str(), 200, 0., 3. );

	
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

      vector<TLorentzVector> reactionVector(reaction->particleList().size());
      Kinematics* kin;
      beamEnergy = cobrem_vs_E->GetRandom();
      ftGen.setBeamEnergy( beamEnergy );  // Resets value of beam energy
      //This while loop will change mass value of subBW specified in arguments
      for(unsigned int iter = 0 ; iter < valVertex.size(); iter++ ){
        // Add BW to associated particle

	// this is the fraction of the central BW distribution that
        // will be generated... throwing away 1% in the tails of
        // the distribution avoids extreme values that cause problems
        // with energy/momentum conservation
        double genFraction = 0.99;
        
        switch( indicateBW[iter] ){
	  case 1: 
            lvMasses[valVertex[iter]] = mpBW[valPDG[iter]](genFraction).first; // Resets value of particle in lower vertex
	    ftGen.setLowerVtxMasses(lvMasses);
	    break;
	  case 2:
            uvMasses[valVertex[iter]] = mpBW[valPDG[iter]](genFraction).first; // Resets value of particle in upper vertex
            ftGen.setUpperVtxMasses(uvMasses);
            break;
        }
      }  
      kin = ftGen.generate(); 
      // Rearranging indices in kinematics class to mimic reactionList
      // Starting with beam
      for(unsigned int k = 0; k < reaction->particleList().size(); k++){
        if( k == 0 ) reactionVector[k] = kin->particle( k );
        if( k > 0 && k <= lvIndices.size()){
	  reactionVector[lvIndices[k-1]] = kin->particle( k ) ;
	}
	if( k > lvIndices.size() && k <= lvIndices.size() + uvIndices.size() ){
	  reactionVector[uvIndices[k - lvIndices.size() - 1]] = kin->particle( k ) ;
	}
      }
      Kinematics* sim = new Kinematics(reactionVector,kin->weight());
      ati.loadEvent( sim, i, batchSize );
      delete kin;
      delete sim;
      i++;
    }
    		
    cout << "Processing events..." << endl;
		
    // include factor of 1.5 to be safe in case we miss peak -- avoid
    // intensity calculation of we are generating flat data
    double maxInten = ( fixedGen ? 0 : 1.5 * ati.processEvents( reaction->reactionName() ) );
		
		
    for( int i = 0; i < batchSize; ++i ){
			
      Kinematics* evt = ati.kinematics( i );
      TLorentzVector recoil;
      for (unsigned int j=0; j < lvIndices.size(); j++){
	      recoil += evt->particle( lvIndices[j] );
       
      }
      TLorentzVector resonance;
      for (unsigned int j= 0; j < uvIndices.size(); j++){
	      resonance += evt->particle( uvIndices[j] );

      }

      // calculate angular variables
      TLorentzVector beam = evt->particle ( 0 );
      TLorentzVector target(0,0,0,targetMass);
      double genWeight = evt->weight();
			
      // cannot ask for the intensity if we haven't called process events above
      double weightedInten = ( fixedGen ? 1 : ati.intensity( i ) ); 
      // cout << " i=" << i << "  intensity_i=" << weightedInten << endl;

      if( !diag ){
				
	// obtain this by looking at the maximum value of intensity * genWeight
	double rand = gRandom->Uniform() * maxInten;
				
	if( weightedInten > rand){


	  m_recoil->Fill( recoil.M() );
	  
	  mW_recoil->Fill( recoil.M(), genWeight );
									
	  t->Fill( -1*(recoil-target).M2() );

   	  tW->Fill( -1*(recoil-target).M2(),genWeight );	  
  
	  e->Fill( beam.E() );

	  eW->Fill( beam.E(),genWeight );

	  intenW->Fill( weightedInten );

	  intenWVsE->Fill( beam.E(), weightedInten );
	  
	  intenWVsM->Fill( resonance.M(), weightedInten );
	
	  m_Meson->Fill( resonance.M() );
	  
	  mW_Meson->Fill( resonance.M(), genWeight );
 
	
	  
	  // we want to save events with weight 1
	  if( !fixedGen ) evt->setWeight( 1.0 );
					
	  if( hddmOut ) hddmOut->writeEvent( *evt, pTypes);
	  rootOut->writeEvent( *evt );
	  ++eventCounter;
	  if(eventCounter >= nEvents) break;
	}
      }
      else{
				
	m_recoil->Fill( recoil.M() );
	mW_recoil->Fill( recoil.M(), genWeight );
	m_Meson->Fill( resonance.M() );
	mW_Meson->Fill( resonance.M(), genWeight );
				
	t->Fill( -1*(recoil-target).M2() );
	tW->Fill( -1*(recoil-target).M2(),genWeight );
	e->Fill( beam.E() );	
	eW->Fill( beam.E(),genWeight );
	eWI->Fill( beam.E(), weightedInten );
	intenW->Fill( weightedInten );
	intenWVsE->Fill( beam.E(), weightedInten );	
	intenWVsM->Fill( resonance.M(), weightedInten );
			
	++eventCounter;
      }
			
      delete evt;
    }
		
    cout << eventCounter << " events were processed." << endl;
  }
	
  m_recoil->Write();
  mW_recoil->Write();
  m_Meson->Write();
  mW_Meson->Write();
  t->Write();
  tW->Write();
  e->Write();
  eW->Write();
  eWI->Write();
  intenW->Write();
  intenWVsE->Write();
  intenWVsM->Write();

  diagOut->Close();
	
  if( hddmOut ) delete hddmOut;
  delete rootOut;
	
  return 0;
}// END OF MAIN()

vector<int> parseString(string vertexString){
  vector<int> Index;
  for(unsigned int i = 0; i < vertexString.size() ; i++){
    Index.push_back(atoi(vertexString.substr(i,1).c_str()));
  }
  return Index;
}// END OF parseString

string checkParticle(string particleString){
  size_t found;
  found = particleString.find( "%" );
  if( found != std::string::npos ){
    particleString = particleString.substr(0, particleString.find( "%" ) );
    return particleString;
  }
  else{
    return particleString;
  }
}//END of checkParticle


inline static unsigned short int UseResonance(Particle_t p)
{
   p = RemapParticleID(p);

	if(IsFixedMass(p) == 1)
		return 0;
	if(p == Unknown)
		return 1;
	if(p == phiMeson)
		return 1;
	if(p == omega)
		return 1;
	return 1;
}//END OF UseResonance
