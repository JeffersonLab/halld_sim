
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "TSystem.h"

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
#include "AMPTOOLS_AMPS/TwoPitdist.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
#include "AMPTOOLS_AMPS/TwoPiAnglesRadiative.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/BreitWigner3body.h"
#include "AMPTOOLS_AMPS/b1piAngAmp.h"
#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/polCoef.h"
#include "AMPTOOLS_AMPS/dblRegge.h"
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
  int numRand = 0;
  
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
    if (arg == "-r"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  numRand = atoi(argv[++i]); }
    if (arg == "-n") useMinos = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
      cout << "   -c <file>\t\t\t\t config file" << endl;
      cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
      cout << "   -r <num random>\t\t\t randomize starting parameters (optional)" << endl;
      exit(1);}
  }
  
  if (configfile.size() == 0){
    cout << "No config file specified" << endl;
    exit(1);
  }

  ConfigFileParser parser(configfile);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();
  string fitName = cfgInfo->fitName();

  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( BreitWigner3body() );
  AmpToolsInterface::registerAmplitude( TwoPSAngles() );
  AmpToolsInterface::registerAmplitude( TwoPSHelicity() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_amp() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_sigma() );
  AmpToolsInterface::registerAmplitude( TwoPitdist() );
  AmpToolsInterface::registerAmplitude( ThreePiAngles() );
  AmpToolsInterface::registerAmplitude( ThreePiAnglesSchilling() );
  AmpToolsInterface::registerAmplitude( TwoPiAnglesRadiative() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( b1piAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( polCoef() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( dblRegge() );
  AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
  
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() ); 
 
  AmpToolsInterface ati( cfgInfo );
  
  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;
  
  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(25000); 
  //fitManager->setPrecision(1e-16);

  // Randomized starting values for multiple fits 
  if(numRand > 0) { 
    
    // keep track of best fit (mininum log-likelihood)
    double minLL = 0;
    int minFitTag = -1;

    // loop over randomized starting values
    for(int ifit=1; ifit <= numRand; ifit++) {
      cout << endl << "###############################" << endl;
      cout << "FIT " << ifit << " OF " << numRand << endl;
      cout << endl << "###############################" << endl;

      ati.randomizeProductionPars(0.5);
      ati.randomizeParameter("dsratio", 0.2, 0.34);

      if( useMinos ) fitManager->minosMinimization();
      else  fitManager->migradMinimization();

      cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

      bool fitFailed = ( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 );
      if( fitFailed ){
	 cout << "ERROR: fit failed use results with caution..." << endl;
      }
    
      // save individual fits
      ati.finalizeFit( Form("%d",ifit) );

      // update best fit
      if( !fitFailed && ati.likelihood() < minLL ) {
	minLL = ati.likelihood();
	minFitTag = ifit;
	//gSystem->Exec(Form("mv %s.ni %s_%d.ni", fitName.data(), fitName.data(), ifit));
      }
    }

    // print best fit results
    if(minFitTag < 0) cout << "ALL FITS FAILED!" << endl;
    else {
      cout << "MINIMUM LIKELIHOOD FROM " << minFitTag << " of " << numRand << " RANDOM PRODUCTION PARS = " << minLL << endl;

      // make best fit available for plotter
      //gSystem->Exec(Form("cp %s_%d.ni %s.ni", fitName.data(), minFitTag, fitName.data()));
      gSystem->Exec(Form("cp %s_%d.fit %s.fit", fitName.data(), minFitTag, fitName.data())); 
    }
  }
  else { // single fit with input parameters from config file

    if( useMinos ) fitManager->minosMinimization();
    else  fitManager->migradMinimization();
    
    bool fitFailed = ( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 );
  
    if( fitFailed ){
      cout << "ERROR: fit failed use results with caution..." << endl;
    }
  
    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;
  
    ati.finalizeFit();
  
    if( seedfile.size() != 0 && !fitFailed ){
      ati.fitResults()->writeSeed( seedfile );
    }
  }
  
  return 0;
}


