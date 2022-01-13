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
#include "AMPTOOLS_AMPS/TwoPiW_brokenetas.h"
#include "AMPTOOLS_AMPS/TwoPitdist.h"
#include "AMPTOOLS_AMPS/TwoPiNC_tdist.h"
#include "AMPTOOLS_AMPS/TwoPiEtas_tdist.h"
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
#include "AMPTOOLS_AMPS/dblReggeMod.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/Piecewise.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, int maxIter, string seedfile) {
   AmpToolsInterface ati( cfgInfo );

   cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

   MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
   fitManager->setMaxIterations(maxIter);

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
      return 1e6;
   }

   cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

   ati.finalizeFit();

   if( seedfile.size() != 0 && !fitFailed ){
      ati.fitResults()->writeSeed( seedfile );
   }

   return ati.likelihood();
}

void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, int maxIter, string seedfile, int numRnd, double maxFraction) {
   AmpToolsInterface ati( cfgInfo );
   string fitName = cfgInfo->fitName();

   cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

   MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
   fitManager->setMaxIterations(maxIter);

   vector< vector<string> > parRangeKeywords = cfgInfo->userKeywordArguments("parRange");

   // keep track of best fit (mininum log-likelihood)
   double minLL = 0;
   int minFitTag = -1;

   for(int i=0; i<numRnd; i++) {
      cout << endl << "###############################" << endl;
      cout << "FIT " << i << " OF " << numRnd << endl;
      cout << endl << "###############################" << endl;

      // randomize parameters
      ati.randomizeProductionPars(maxFraction);
      for(size_t ipar=0; ipar<parRangeKeywords.size(); ipar++) {
         ati.randomizeParameter(parRangeKeywords[ipar][0], atof(parRangeKeywords[ipar][1].c_str()), atof(parRangeKeywords[ipar][2].c_str()));
      }

      if(useMinos)
         fitManager->minosMinimization();
      else
         fitManager->migradMinimization();

      bool fitFailed = (fitManager->status() != 0 && fitManager->eMatrixStatus() != 3);

      if( fitFailed )
         cout << "ERROR: fit failed use results with caution..." << endl;

      cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

      ati.finalizeFit(to_string(i));

      if( seedfile.size() != 0 && !fitFailed ){
         string seedfile_rand = seedfile + Form("_%d.txt", i);
         ati.fitResults()->writeSeed( seedfile_rand );
      }

      // update best fit
      if( !fitFailed && ati.likelihood() < minLL ) {
         minLL = ati.likelihood();
         minFitTag = i;
      }
   }

   // print best fit results
   if(minFitTag < 0) cout << "ALL FITS FAILED!" << endl;
   else {
      cout << "MINIMUM LIKELIHOOD FROM " << minFitTag << " of " << numRnd << " RANDOM PRODUCTION PARS = " << minLL << endl;
      gSystem->Exec(Form("cp %s_%d.fit %s.fit", fitName.data(), minFitTag, fitName.data()));
      if( seedfile.size() != 0 )
         gSystem->Exec(Form("cp %s_%d.txt %s.txt", seedfile.data(), minFitTag, seedfile.data()));
   }
}

void runParScan(ConfigurationInfo* cfgInfo, bool useMinos, int maxIter, string seedfile, string parScan) {
   double minVal, maxVal, stepSize;
   int steps;

   vector< vector<string> > parScanKeywords = cfgInfo->userKeywordArguments("parScan");

   if(parScanKeywords.size()==0) {
      cout << "No parScan keyword found in configuration file. Set up at least one parameter for scanning! Aborting." << endl;
      return;
   } else {
      for(size_t ipar=0; ipar<parScanKeywords.size(); ipar++) {
         if(parScanKeywords[ipar][0]==parScan) {
            minVal = atof(parScanKeywords[ipar][1].c_str());
            maxVal = atof(parScanKeywords[ipar][2].c_str());
            stepSize = atof(parScanKeywords[ipar][3].c_str());
            steps = trunc((maxVal-minVal)/stepSize)+1;
            break;
         } else
            cout << "Skipping configuration to scan " << parScanKeywords[ipar][0] << "since scanning of " << parScan << " was requested..." << endl;
      }
   }

   AmpToolsInterface ati( cfgInfo );

   string fitName = cfgInfo->fitName();
   cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

   ParameterManager* parMgr = ati.parameterManager();
   MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
   fitManager->setMaxIterations(maxIter);


   for(int i=0; i<steps; i++) {
      cout << endl << "###############################" << endl;
      cout << "FIT " << i << " OF " << steps << endl;
      cout << endl << "###############################" << endl;

      // set parameter to be scanned
      vector<ParameterInfo*> parInfoVec = cfgInfo->parameterList();

      auto parItr = parInfoVec.begin();
      for( ; parItr != parInfoVec.end(); ++parItr ) {
         if( (**parItr).parName() == parScan ) break;
      }

      if( parItr == parInfoVec.end() ){
         cout << "ERROR:  request to scan nonexistent parameter:  " << parScan << endl;
         return;
      }

      // set and fix parameter for scan
      double value = minVal + i*stepSize;
      parMgr->setAmpParameter( parScan, value );

      cfgInfo->setFitName(fitName + "_scan");

      if(useMinos)
         fitManager->minosMinimization();
      else
         fitManager->migradMinimization();

      bool fitFailed = (fitManager->status() != 0 && fitManager->eMatrixStatus() != 3);

      if( fitFailed )
         cout << "ERROR: fit failed use results with caution..." << endl;

      cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

      ati.finalizeFit(to_string(i));

      if( seedfile.size() != 0 && !fitFailed ){
         string seedfile_scan = seedfile + Form("_scan_%d.txt", i);
         ati.fitResults()->writeSeed( seedfile_scan );
      }
   }
}

int main( int argc, char* argv[] ){

   // set default parameters

   bool useMinos = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   int maxIter = 10000;

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
         else  numRnd = atoi(argv[++i]); }
      if (arg == "-m"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  maxIter = atoi(argv[++i]); }
      if (arg == "-n") useMinos = true;
      if (arg == "-p"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  scanPar = argv[++i]; }
      if (arg == "-h"){
         cout << endl << " Usage for: " << argv[0] << endl << endl;
         cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
         cout << "   -c <file>\t\t\t\t config file" << endl;
         cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
         cout << "   -r <int>\t\t\t Perform <int> fits each seeded with random parameters" << endl;
         cout << "   -p <parameter> \t\t\t\t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
         cout << "   -m <int>\t\t\t Maximum number of fit iterations" << endl; 
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
   AmpToolsInterface::registerAmplitude( Piecewise() );

   AmpToolsInterface::registerDataReader( ROOTDataReader() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );

   if(numRnd==0){
      if(scanPar=="")
         runSingleFit(cfgInfo, useMinos, maxIter, seedfile);
      else
         runParScan(cfgInfo, useMinos, maxIter, seedfile, scanPar);
   } else {
      runRndFits(cfgInfo, useMinos, maxIter, seedfile, numRnd, 0.5);
   }

   return 0;
}


