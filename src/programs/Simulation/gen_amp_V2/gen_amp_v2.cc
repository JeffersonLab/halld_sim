
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
int main( int argc, char* argv[] ){
  
  TString beamConfigFile("");
  string  atConfigFile("");
  string  hddmName("");
  string  outName("");
  string  lvString("");
  string  uvString("");

	
  bool genFlat = false;
  bool fsRootFormat = false;
  bool diag = false;
  bool uvPSRange = false;
  bool lvPSRange = false;
  bool tRange = false;
  vector<double> uvBW;
  vector<double> lvBW;
  vector<double> tSlope;

  // default upper and lower bounds -- these
  // just need to be outside the kinematic bounds
  // of any reaction
  double lvMin = 0.0;
  double uvMin = 0.0;
  double lvMax = 5.0;
  double uvMax = 5.0;
  
  double beamMaxE   = 12.0;
  double beamPeakE  = 9.0;
  double beamLowE   = 3.0;
  double beamHighE  = 12.0;

  int runNum = 30731;
  unsigned int seed = 0;
  unsigned int reWeight =7;
  double tMin = 0.0;
  double tMax = 15.0;

  int nEvents = 10000;
  int batchSize = 10000;
  int uvCounter = 0;
  int lvCounter = 0;
  int tSlopeCounter = 0;

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
	uvBW.push_back( atof( argv[i+1] ) ); 
        uvBW.push_back( atof( argv[i+2] ) ); 
        uvBW.push_back( atof( argv[i+3] ) ); 
	uvCounter++;
	
      }
    }
    if (arg == "-lvBW") {
      if ((i+1 == argc) || (i+2 == argc) || (i+3 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-') || (argv[i+3][0] == '-')) arg = "-h";
      else{ 
        lvBW.push_back( atof( argv[i+1] ) ); 
        lvBW.push_back( atof( argv[i+2] ) ); 
        lvBW.push_back( atof( argv[i+3] ) ); 
	lvCounter++;
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
	lvMin = atof( argv[i+1] ); 
        lvMax = atof( argv[i+2] ); 
        lvPSRange = true;
      }
    }
    if (arg == "-uvRange"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
	uvMin = atof( argv[i+1] ); 
        uvMax = atof( argv[i+2] ); 
        uvPSRange = true;
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
	tSlope.push_back( atof( argv[i+1] ));
	tSlope.push_back( atof( argv[i+2] ));
	tSlopeCounter++;
      }
    }
    if (arg == "-tRange"){
      if ((i+1 == argc) || (i+2 == argc) || (argv[i+1][0] == '-') || (argv[i+2][0] == '-')) arg = "-h";
      else{  
	tMin = atof( argv[i+1] ); 
        tMax = atof( argv[i+2] ); 
        tRange = true;
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
    if (arg == "-f") genFlat = true; 
    if (arg == "-d") diag = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -ac      <file>\t AmpTools config file" << endl;
      cout << "\t -bc      <file>\t beam config file (will generate local file within coherent peak if not included) " << endl;
      cout << "\t -o       <name>\t ROOT file output name" << endl;
      cout << "\t -uv      <string>\t upper-vertex indices (0 is first) in AmpTools reaction" << endl;
      cout << "\t -lv      <string>\t lower-vertex indices (0 is first) in AmpTools reaction" << endl;
      cout << "\t -fsroot          \t Enable output in FSRoot format [optional]" << endl;
      cout << "\t -hd      <name>\t HDDM file output name [optional]" << endl;
      cout << "\t -f              \t Generate flat in phase space [optional]" << endl;
      cout << "\t -uvBW    <value>\t concentrate events around one or more upper vertex resonances (Mass Width Scale) [optional]" << endl; 
      cout << "\t -lvBW    <value>\t concentrate events around one or more lower vertex resonances (Mass Width Scale) [optional]" << endl; 
      cout << "\t -n       <value>\t Minimum number of events to generate [optional]" << endl;
      cout << "\t -r       <value>\t Run number assigned to generated events [optional]" << endl;
      cout << "\t -s       <value>\t Random number seed initialization [optional]" << endl;
      cout << "\t -lvRange <value>\t Lower vertex mass range ( lowerLimit upperLimit ) (GeV) [optional]" << endl;
      cout << "\t -uvRange <value>\t Upper vertex mass range  ( lowerLimit upperLimit ) (GeV) [optional]" << endl;
      cout << "\t -m       <value>\t Electron beam energy (or photon energy endpoint) [optional]" << endl;
      cout << "\t -p       <value>\t Coherent peak photon energy [optional]" << endl;
      cout << "\t -a       <value>\t Minimum photon energy to simulate events [optional]" << endl;
      cout << "\t -b       <value>\t Maximum photon energy to simulate events [optional]" << endl;
      cout << "\t -t       <value>\t Momentum transfer slope ( t Scale ) [optional]" << endl;
      cout << "\t -tRange  <value>\t momentum transfer range (lowerLimit upperLimit) (Gev)^2 [optional]" << endl;
      cout << "\t -d              \t Plot only diagnostic histograms [optional]" << endl; 
      cout << "\t -mask    <value>\t Remove weight to recover phase space for upper vertex mass, lower vertex, or momentum transfer t respectively ( uvM lvM t ) [optional]" << endl << endl; 
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
  /// get the masses for each.....

  vector< double > lvMasses;
  vector< double > uvMasses;
  vector< int > pTypes;
  ostringstream locStream;
  ostringstream locRecoilStream;
 
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
	  locRecoilStream << ParticleName_ROOT( particle );
        }
      }
      // this loop will check if particle is part of upper vertex according to user
      for( unsigned int k = 0; k < uvIndices.size(); ++k){
        if( uvIndices[k] > 0 && ((unsigned int)uvIndices[k] == i )){
          cout << "This is upper vertex particle with indices " << uvIndices[k] 
          << " with name " <<  tempString << endl;
          uvMasses.push_back( ParticleMass( particle ) ); 
          locStream << ParticleName_ROOT( particle );     
        }
      }
      // add particle to pTypes
      pTypes.push_back( particle );
    }
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

  FixedTargetGenerator ftGen( beamEnergy, targetMass, uvMasses, lvMasses );

  
  // Include optional kinematic ranges (working on including BW values)
  if( lvPSRange ){ 
    ftGen.setUpperVtxRange( uvMin, uvMax ); }
  if( uvPSRange ){
    ftGen.setLowerVtxRange( lvMin, lvMax ); }
  if( tRange ){
    ftGen.setMomentumTransferRange( tMin, tMax ); }
  for(int i=0; i<tSlopeCounter; i++){
    ftGen.addMomentumTransfer( tSlope[2*i], tSlope[2*i + 1] ); 
  }
   
  // the lines below should be tailored by the user for the particular desired
  // set of amplitudes -- doing so will improve efficiency.  Leaving as is
  // won't make MC incorrect, it just won't be as fast as it could be
		
  for(int i=0; i< uvCounter; i++){
    ftGen.addUpperVtxBW( uvBW[3*i] , uvBW[3*i + 1] , uvBW[3*i + 2] );
  }
  for(int i=0; i< lvCounter; i++){ 
    ftGen.addLowerVtxBW( lvBW[3*i] , lvBW[3*i + 1] , lvBW[3*i + 2] );
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
  
  string locHistTitle = string("Meson Mass ;") + locStream.str() + string(" Invariant Mass (GeV/c^{2});");
  string locRecoilTitle = string("Baryon Mass ;") + locRecoilStream.str() + string(" Invariant Mass (GeV/c^{2});");

	
  TH1F* t = new TH1F( "t", "-t Distribution", 200, 0, 2 );
  TH1F* tW = new TH1F( "tW", "-t Distribution", 200, 0, 2 ); 
  TH1F* E = new TH1F( "E", "Beam Energy", 120, 0, 12 );
  TH1F* EW = new TH1F( "EW", "Beam Energy", 120, 0, 12 );
  TH1F* EWI = new TH1F( "EWI", "Beam Energy", 120, 0, 12 );
  TH1F* intenW = new TH1F("intenW", "True PDF/ Gen. PDF", 1000, -0.1e-03, 0.8e-03);
  TH2F* intenWVsE = new TH2F("intenWVsE","Ratio vs. M", 100, 0, 12, 1000, -0.1e-03, 0.8e-03);
  
  TH1F* M_Meson = new TH1F( "M_Meson", locHistTitle.c_str(), 200, uvMin, uvMax );
  TH1F* MW_Meson = new TH1F( "MW_Meson", locHistTitle.c_str(), 200, uvMin, uvMax );
  TH1F* M_recoil = new TH1F( "M_recoil", locRecoilTitle.c_str(), 200, lvMin, lvMax );
  TH1F* MW_recoil = new TH1F( "MW_recoil", locRecoilTitle.c_str(), 200, lvMin, lvMax );

	
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
      ftGen.setBeamEnergy( cobrem_vs_E->GetRandom() );  // Resets value of beam energy
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
    double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
		
		
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
      double weightedInten = ( genFlat ? 1 : ati.intensity( i ) ); 
      // cout << " i=" << i << "  intensity_i=" << weightedInten << endl;

      if( !diag ){
				
	// obtain this by looking at the maximum value of intensity * genWeight
	double rand = gRandom->Uniform() * maxInten;
				
	if( weightedInten > rand || genFlat ){


	  M_recoil->Fill( recoil.M() );
	  
	  MW_recoil->Fill( recoil.M(), genWeight );
									
	  t->Fill( -1*(recoil-target).M2() );

   	  tW->Fill( -1*(recoil-target).M2(),genWeight );	  
  
	  E->Fill( beam.E() );

	  EW->Fill( beam.E(),genWeight );

	  intenW->Fill( weightedInten );

	  intenWVsE->Fill( beam.E(), weightedInten );
	
	  M_Meson->Fill( resonance.M() );
	  
	  MW_Meson->Fill( resonance.M(), genWeight );
 
	
	  
	  // we want to save events with weight 1
	  evt->setWeight( 1.0 );
					
	  if( hddmOut ) hddmOut->writeEvent( *evt, pTypes);
	  rootOut->writeEvent( *evt );
	  ++eventCounter;
	  if(eventCounter >= nEvents) break;
	}
      }
      else{
				
	M_recoil->Fill( recoil.M() );
	MW_recoil->Fill( recoil.M(), genWeight );
	M_Meson->Fill( resonance.M() );
	MW_Meson->Fill( resonance.M(), genWeight );
				
	t->Fill( -1*(recoil-target).M2() );
	tW->Fill( -1*(recoil-target).M2(),genWeight );
	E->Fill( beam.E() );	
	EW->Fill( beam.E(),genWeight );
	EWI->Fill( beam.E(), weightedInten );
	intenW->Fill( weightedInten );
	intenWVsE->Fill( beam.E(), weightedInten );	
			
	++eventCounter;
      }
			
      delete evt;
    }
		
    cout << eventCounter << " events were processed." << endl;
  }
	
  M_recoil->Write();
  MW_recoil->Write();
  M_Meson->Write();
  MW_Meson->Write();
  t->Write();
  tW->Write();
  E->Write();
  EW->Write();
  EWI->Write();
  intenW->Write();
  intenWVsE->Write();

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
