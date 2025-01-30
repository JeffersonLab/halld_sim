
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

#include "GLUEX_EVTGEN/EvtGenDecayer.h"

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


bool compVector(int a, int b);
pair< vector<int>,vector<bool> > parseString(string vertexString);
vector<string> trueReactionDecay( vector<string>& keywordArgs, vector<string> reactionList, vector<int> uvIndices, vector<int> lvIndices, pair< vector<int>,vector<bool> >& orderTrueuvIndices, pair< vector<int>,vector<bool> >& orderTruelvIndices);
vector<TLorentzVector> makeReaction( Kinematics* kin, vector< string > keywordArgs, vector< TLorentzVector > reactionVector, vector< int > uvIndex, vector< int > lvIndex, int numUvDecay, int numLvDecay );
string checkParticle(string particleString);
vector<double> getVertexMasses(vector<string>& reactionList, map<string, BreitWignerGenerator>& mpBWGen, vector<int>& indices, vector<bool>& hasResonance);
vector<int> getTypes(vector<string>& reactionList);

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
  
  map<string, BreitWignerGenerator> mpBW;
  // the map index must match what is written in the AmpTools config file
  // the particle enum lookup needs to match what is in particleType.h
  mpBW["Omega"] = BreitWignerGenerator( ParticleMass( (Particle_t)omega ), 0.00868, seed); // Initialize BW for omega
  mpBW["Phi"] = BreitWignerGenerator( ParticleMass( (Particle_t)phiMeson ), 0.004249, seed); // Initialize BW for omega

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
  vector< string > pList = reaction->particleList();
  pair< vector< int >,vector<bool> > lvIndices = parseString( lvString );
  pair< vector< int >,vector<bool> > uvIndices = parseString( uvString );
  
  // This section initiates variables needed for EvtGen decay 
  vector< string > truepList;
  pair< vector< int >,vector<bool> > orderTruelvIndices;
  pair< vector< int >,vector<bool> > orderTrueuvIndices;
  EvtGenDecayer* decayer;
  vector< vector<string> > trueReactionKeyword = cfgInfo->userKeywordArguments("trueReaction");
  vector< string > keywordArgs = trueReactionKeyword[0];
  truepList = reaction->particleList(); 
  

  // Redefinition of variables if keyword is found
  if(trueReactionKeyword.size() == 1){
  truepList = trueReactionDecay( keywordArgs, pList, uvIndices.first, lvIndices.first, orderTrueuvIndices, orderTruelvIndices);
  }

  vector< int > pTypes;
  ostringstream upperVtxName;
  ostringstream lowerVtxName;

  vector< double > uvMasses;
  vector< double > lvMasses;
  uvMasses = getVertexMasses( pList, mpBW, uvIndices.first, uvIndices.second );
  lvMasses = getVertexMasses( pList, mpBW, lvIndices.first, lvIndices.second );
  if(trueReactionKeyword.size() == 1){
  uvMasses = getVertexMasses( truepList, mpBW, orderTrueuvIndices.first, orderTrueuvIndices.second );
  lvMasses = getVertexMasses( truepList, mpBW, orderTruelvIndices.first, orderTruelvIndices.second );	  
  }

  pTypes = getTypes( pList );

   // This section will generate daughter particles for specified decay particles
  if( trueReactionKeyword.size() == 1 ){
    // First check local decay config file if found
    ifstream file("userDecay.dec");
    if( !file.good() ){
      cout << "ERROR:  Missing local EvtGen decay config file" << endl;
      exit( 1 );
    }
    decayer = new EvtGenDecayer();
  }

  
  
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
  for( map< string, BreitWignerGenerator >::iterator mapItr = mpBW.begin(); mapItr != mpBW.end(); ++mapItr ){
    mapItr->second.setSeed( seed );
  } 

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
  
  string locHistTitle = string("Meson Mass ; Upper Vertex") + string(" Invariant Mass (GeV/c^{2});");
  string locRecoilTitle = string("Baryon Mass ; Lower Vertex") + string(" Invariant Mass (GeV/c^{2});");

	
  TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );
  TH1F* tW = new TH1F( "tW", "-t Distribution", 200, 0, 2 ); 
  TH1F* e = new TH1F( "e", "Beam Energy", 120, 0, 12 );
  TH1F* eW = new TH1F( "eW", "Beam Energy", 120, 0, 12 );
  TH1F* eWI = new TH1F( "eWI", "Beam Energy", 120, 0, 12 );
  TH1F* intenW = new TH1F("intenW", "True PDF/ Gen. PDF", 1000, -0.1e-03, 0.8e-03);
  TH2F* intenWVsE = new TH2F("intenWVsE","Ratio vs. E", 100, 0, 12, 200, -0.1e-03, 0.8e-03);
  TH2F* intenWVsM = new TH2F("intenWVsM","Ratio vs. M", 100, 0, 3., 200, -0.1e-03, 0.8e-03);
  
  TH1F* m_Meson = new TH1F( "m_Meson", locHistTitle.c_str(), 200, 0., 3. );
  TH1F* mW_Meson = new TH1F( "mW_Meson", locHistTitle.c_str(), 200, 0., 3. );
  TH1F* m_recoil = new TH1F( "m_recoil", locRecoilTitle.c_str(), 200, 0., 3. );
  TH1F* mW_recoil = new TH1F( "mW_recoil", locRecoilTitle.c_str(), 200, 0., 3. );


  // Making temporary uv/lv indices for reorganization
  pair< vector< int >,vector<bool> > tempuvIndices = ( trueReactionKeyword.size() == 1 ? orderTrueuvIndices : uvIndices );	
  pair< vector< int >,vector<bool> > templvIndices = ( trueReactionKeyword.size() == 1 ? orderTruelvIndices : lvIndices );	
  vector< string > temppList = ( trueReactionKeyword.size() == 1 ? truepList : pList);
  
  
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
      //This will change mass value of upper/lower vertex if particle is an omega or phi
      if( count(tempuvIndices.second.begin(), tempuvIndices.second.end(), true) > 0  ){
        ftGen.setUpperVtxMasses( getVertexMasses( temppList, mpBW, tempuvIndices.first, tempuvIndices.second ) );
      }
      if( count(templvIndices.second.begin(), templvIndices.second.end(), true) > 0 ) ftGen.setUpperVtxMasses( getVertexMasses( temppList, mpBW, templvIndices.first, templvIndices.second ) );
      kin = ftGen.generate(); 
    
      
      for(int h = 0; h < kin->particleList().size(); h++){
          kin->particle( h ).Print();
          cout << endl << endl;
      }

      // This section will generate daughter particles for specified decay particles
      if( trueReactionKeyword.size() == 1 ){
  	int uvIndex = 0; // Quick bookeeping of upper index
        int lvIndex = 0; // Qucik bookeeping of lower index
	vector<TLorentzVector> trueReactionVector;
	trueReactionVector.push_back( kin->particle(0) ); //Add Beam 4-vector to new vector
	for( int h = 0; h < (int)keywordArgs.size()/3; h++){
	  if( keywordArgs[3*h] == "lv" ){
            cout << " lv particle to be decayed " << checkParticle( temppList[ templvIndices.first[lvIndex] ] ).c_str() <<  " for event " << i << endl;
	    vector< pair<TLorentzVector, int> > children = decayer->decayParticle( kin->particle( h + 1 ), ParticleEnum( checkParticle( temppList[ templvIndices.first[lvIndex] ] ).c_str() ) );
            lvIndex++;
	    for( auto child_itr = children.begin(); child_itr != children.end(); child_itr++){
            trueReactionVector.push_back( child_itr->first );
	    }
          } 
		
	  if( keywordArgs[3*h] == "uv" ){ 
  	    cout << "uv particle to be decayed " << checkParticle( temppList[ tempuvIndices.first[uvIndex] ] ).c_str() << endl << " for event " << i << endl;
	    vector< pair<TLorentzVector, int> > children = decayer->decayParticle( kin->particle( templvIndices.first.size() + h + 1 ), ParticleEnum( checkParticle( temppList[ tempuvIndices.first[uvIndex] ] ).c_str() ) );
	    uvIndex++;
	    for( auto child_itr = children.begin(); child_itr != children.end(); child_itr++){
            trueReactionVector.push_back( child_itr->first );
            }
	  }
        }
	for( int g = 0; g < trueReactionVector.size(); g++){
          trueReactionVector[g].Print();
	}
	cout << "Particles coming out of EvtGen is " << trueReactionVector.size() << endl;
	vector<TLorentzVector> trialVector = makeReaction( kin, temppList, trueReactionVector, tempuvIndices.first, templvIndices.first, uvIndex, lvIndex);
	vector<TLorentzVector> trialVector2 = makeReaction( kin, temppList, trueReactionVector, tempuvIndices.first, templvIndices.first, uvIndex, lvIndex);
        kin->setParticleList( trialVector2 );
        //kin->setParticleList( makeReaction( kin, temppList, trueReactionVector, tempuvIndices.first, templvIndices.first, uvIndex, lvIndex) );
	for(int g = 0; g < kin->particleList().size(); g++){
          kin->particle( g ).Print();
          trialVector[g].Print();
          trialVector2[g].Print();
	  cout << endl << endl;
        }
	cout << "Passed the stuff " << endl;
      }  

      // Rearranging indices in kinematics class to mimic reactionList
      // Starting with beam
      for(unsigned int k = 0; k < reaction->particleList().size(); k++){
        if( k == 0 ){ //This is for the beam
          reactionVector[k] = kin->particle( k );
        }
        if( k > 0 && k <= lvIndices.first.size()){ //This is for the lower vertex 
	  reactionVector[lvIndices.first[k-1]] = kin->particle( k ) ;
	}
	if( k > lvIndices.first.size() && k <= lvIndices.first.size() + uvIndices.first.size() ){//This is for the upper vertex
	  reactionVector[uvIndices.first[k - lvIndices.first.size() - 1]] = kin->particle( k ) ;
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
      for (unsigned int j=0; j < lvIndices.first.size(); j++){
	      recoil += evt->particle( lvIndices.first[j] );
       
      }
      TLorentzVector resonance;
      for (unsigned int j= 0; j < uvIndices.first.size(); j++){
	      resonance += evt->particle( uvIndices.first[j] );

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
  if( trueReactionKeyword.size() == 1 ) delete decayer;

  return 0;
}// END OF MAIN()

// Comparator function
bool compVector(int a, int b){
  return a < b;

}//END of compVector

pair< vector<int>,vector<bool> > parseString(string vertexString){
  vector<int> index;
  vector<bool> boolDex;
  for(unsigned int i = 0; i < vertexString.size() ; i++){
    index.push_back( atoi(vertexString.substr(i,1).c_str()) );
    boolDex.push_back( false );
  }
  sort(index.begin(), index.end(), compVector);
  return make_pair( index,boolDex );
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


// This function will create a new reactionList order and names if user specifies a decay using the Keyword "trueReaction"
// First step is to replace beginning indice of the reactionList with the parent particle specified
// Then remove the remaining daughter particles from reactionList, this will be passed through FixedTargetGenerator
// Lastly, we need to generate a new index vector for the upper/lower vertexes  

vector<string> trueReactionDecay(vector<string>& keywordArgs, vector<string> reactionList, vector<int> uvIndices, vector<int> lvIndices, pair< vector<int>,vector<bool> >& orderTrueuvIndices, pair< vector<int>,vector<bool> >& orderTruelvIndices){

vector<int> tempOrderTrueuvIndicesFirst;
  vector<bool> tempOrderTrueuvIndicesSecond;
  vector<int> tempOrderTruelvIndicesFirst;
  vector<bool> tempOrderTruelvIndicesSecond;
  vector<string> tempReactionList;
  vector<int> tempuvIndex = uvIndices;
  vector<int> templvIndex = lvIndices;
  int tempInt = 0;
  //Need the next two ints to keep track of total number of true uv/lv particles
  int tempuvInt = 0;
  int templvInt = 0;

  tempReactionList.push_back(reactionList[0]); // Add beam index to list 
  for( int k = 0; k < (int)keywordArgs.size()/3; k++){
    tempInt++;
    pair< vector<int>, vector<bool> > trueIndices = parseString( keywordArgs[3*k+2] );
    // First add parent name to tempReactionList        
    tempReactionList.push_back(keywordArgs[3*k+1]);
    // Add particle to index vector 
    if( keywordArgs[3*k] == "uv" ){
      tempOrderTrueuvIndicesFirst.push_back( tempInt );
      tempOrderTrueuvIndicesSecond.push_back( false );
    }
    if( keywordArgs[3*k] == "lv" ){
      tempOrderTruelvIndicesFirst.push_back( tempInt );
      tempOrderTruelvIndicesSecond.push_back( false );
    }
  }
  // Remove decay indices from uv/lv array. Will be used for non-decaying particles
  for( int k = 0; k < (int)keywordArgs.size()/3; k++){
    pair< vector<int>, vector<bool> > trueIndices = parseString( keywordArgs[3*k+2] );
    for( int j = 0; j < (int)trueIndices.first.size(); j++){
      if( keywordArgs[3*k] == "uv" ){
        tempuvIndex.erase( find( tempuvIndex.begin(), tempuvIndex.end(), trueIndices.first[j] ) );
        tempuvInt++;
      }
      if( keywordArgs[3*k] == "lv" ){
        templvIndex.erase( find( templvIndex.begin(), templvIndex.end(), trueIndices.first[j] ) );
        templvInt++;
      }
    }
  }
  tempInt++; //This is to make sure index is correct
  // Add the non-decaying particles into new uv/lv array
  for( int m = 0; m < (int)uvIndices.size() - tempuvInt; m++){
    tempOrderTrueuvIndicesFirst.push_back( tempInt );
    tempOrderTrueuvIndicesSecond.push_back( false );
    tempInt++;
    tempReactionList.push_back( reactionList[ tempuvIndex[m] ] );
  }

  for( int v = 0; v < (int)lvIndices.size() - templvInt; v++){
    tempOrderTruelvIndicesFirst.push_back( tempInt );
    tempOrderTruelvIndicesSecond.push_back( false );
    tempReactionList.push_back( reactionList[ templvIndex[v] ] );
    tempInt++;
  }

  // This is needed to add to NULL pair
  orderTrueuvIndices.first.resize(tempOrderTrueuvIndicesFirst.size());
  orderTrueuvIndices.second.resize(tempOrderTrueuvIndicesSecond.size());
  orderTruelvIndices.first.resize(tempOrderTruelvIndicesFirst.size());
  orderTruelvIndices.second.resize(tempOrderTruelvIndicesSecond.size());


  for( int j = 0; j < (int)tempOrderTrueuvIndicesFirst.size(); j++){
    orderTrueuvIndices.first[j] = tempOrderTrueuvIndicesFirst[j];
    orderTrueuvIndices.second[j] = tempOrderTrueuvIndicesSecond[j];
    cout << orderTrueuvIndices.first[j] << endl;
  }
  for( int m = 0; m < (int)tempOrderTruelvIndicesFirst.size(); m++){
   orderTruelvIndices.first[m] = tempOrderTruelvIndicesFirst[m];
   orderTruelvIndices.second[m] = tempOrderTruelvIndicesSecond[m];
   cout << orderTruelvIndices.first[m] << endl;
  }
  return tempReactionList;

  
}//END of trueReactionDecay 	

vector< TLorentzVector > makeReaction( Kinematics* kin, vector<string> keywordArgs, vector<TLorentzVector> reactionVector, vector<int> uvIndex, vector<int> lvIndex, int numUvDecay, int numLvDecay ){

  vector< TLorentzVector > tempLorentz = reactionVector;
  int lowerVNum = 1; // Statrting with 1 cause the beam is included
  int upperVNum = 0;
  // First thing to do is find where to add the non-decayed particles since the decayed particles are already in tempLorentz
  for( int i = 0; i < (int)keywordArgs.size()/3; i++ ){
    if( keywordArgs[3*i] == "lv" ) lowerVNum += keywordArgs[3*i +2].length();
    if( keywordArgs[3*i] == "uv" ) upperVNum += keywordArgs[3*i +2].length();
  }
  // Add the non-decaying particles to tempLorentz starting with the lower vertex
  for( int j = 0; j < (int)lvIndex.size() - numLvDecay; j++ ){
    tempLorentz.insert( tempLorentz.begin() + lowerVNum + j, kin->particle( numLvDecay + j + 1) );
  }
  // The non-decaying uv particles go at the end
  for( int k = 0; k < (int)uvIndex.size() - numUvDecay; k++ ){
    tempLorentz.push_back( kin->particle( lvIndex.size() + numUvDecay + k + 1) );
  }
  return tempLorentz;

}// END OF makeReaction



vector<double> getVertexMasses(vector<string>& reactionList, map<string, BreitWignerGenerator>& mpBWGen, vector<int>& indices, vector<bool>& hasResonance){
  vector<double> vertexMasses;
  // this is the fraction of the central BW distribution that
  // will be generated... throwing away 1% in the tails of
  // the distribution avoids extreme values that cause problems
  // with energy/momentum conservationdouble 
  double genFraction = 0.99;
  for(unsigned int i = 0; i < indices.size(); i++){
    string tempString = checkParticle( reactionList[indices[i]] );
    Particle_t particle = ParticleEnum( tempString.c_str() ); 
    if (particle == 0 && i > 0){
      cout << "ERROR:  unknown particle " << tempString 
           << " unable to configure generator." << endl;
      exit( 1 );
    }
    else{
      vertexMasses.push_back( mpBWGen.find(tempString) == mpBWGen.end() ? ParticleMass( particle ) : mpBWGen[tempString](genFraction).first );
      hasResonance[i] = mpBWGen.find(tempString) != mpBWGen.end();
    }
  }
  return vertexMasses;
 
}//END OF getVertexMasses

vector<int> getTypes(vector<string>& reactionList){
  vector<int> pTypes;
  for(unsigned int i = 0; i < reactionList.size(); i++){
    string tempString = checkParticle( reactionList[i] );
    Particle_t particle = ParticleEnum( tempString.c_str() );
    if (particle == 0 && i > 0){
      cout << "ERROR:  unknown particle " << tempString
           << " unable to configure generator." << endl;
      exit( 1 );
    }
    else{
      pTypes.push_back(particle); 
    }
  }
  return pTypes;
}//END OF getTypes

