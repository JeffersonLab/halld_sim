
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderWithTCut.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
#include "AMPTOOLS_AMPS/TwoPiWt_primakoff.h"
#include "AMPTOOLS_AMPS/TwoPiWt_sigma.h"
#include "AMPTOOLS_AMPS/TwoPiW_brokenetas.h"
#include "AMPTOOLS_AMPS/TwoPitdist.h"
#include "AMPTOOLS_AMPS/TwoPiNC_tdist.h"
#include "AMPTOOLS_AMPS/TwoPiEtas_tdist.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
#include "AMPTOOLS_AMPS/VecRadiative_SDME.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/BreitWigner3body.h"
#include "AMPTOOLS_AMPS/b1piAngAmp.h"
#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/polCoef.h"
#include "AMPTOOLS_AMPS/dblRegge.h"
#include "AMPTOOLS_AMPS/dblReggeMod.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
	
  // set default parameters
  
  bool useMinos = false;

  string configfile;
  string seedfile;
  
  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-s"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  seedfile = argv[++i]; }
    if (arg == "-n") useMinos = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
      cout << "   -c <file>\t\t\t\t config file" << endl;
      cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
      exit(1);}
  }
  
  if (configfile.size() == 0){
    cout << "No config file specified" << endl;
    exit(1);
  }

  ConfigFileParser parser(configfile);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( BreitWigner3body() );
  AmpToolsInterface::registerAmplitude( TwoPSAngles() );
  AmpToolsInterface::registerAmplitude( TwoPSHelicity() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_amp() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_sigma() );
  AmpToolsInterface::registerAmplitude( TwoPiW_brokenetas() );
  AmpToolsInterface::registerAmplitude( TwoPitdist() );
  AmpToolsInterface::registerAmplitude( TwoPiNC_tdist() );
  AmpToolsInterface::registerAmplitude( TwoPiEtas_tdist() );
  AmpToolsInterface::registerAmplitude( ThreePiAngles() );
  AmpToolsInterface::registerAmplitude( ThreePiAnglesSchilling() );
  AmpToolsInterface::registerAmplitude( VecRadiative_SDME() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( b1piAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( polCoef() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( dblRegge() );
  AmpToolsInterface::registerAmplitude( dblReggeMod() );
  AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
  
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() ); 
 
  AmpToolsInterface ati( cfgInfo );

  double likelihood_init =  ati.likelihood();
  
  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << likelihood_init << endl;
  if (!isfinite(likelihood_init)) {
      cout << "*** fit- infinite Likelihood before minimization EXIT***" << endl;
      exit(1);
    }
  
  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(10000); 
 
  if( useMinos ){
    
    fitManager->minosMinimization();
  }
  else{
    
    fitManager->migradMinimization();
  }
  
  bool fitFailed =
    ( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 );
  
  if( fitFailed ){
    cout << "ERROR: fit failed use results with caution..." << endl;
  }
  
  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;
  
  ati.finalizeFit();
  
  if( seedfile.size() != 0 && !fitFailed ){
    
    ati.fitResults()->writeSeed( seedfile );
  }
  
  return 0;
}


