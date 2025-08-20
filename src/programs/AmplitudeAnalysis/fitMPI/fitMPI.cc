#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <limits>

#include "TSystem.h"

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderWithTCut.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
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
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/polCoef.h"
#include "AMPTOOLS_AMPS/DblRegge_FastEta.h"
#include "AMPTOOLS_AMPS/DblRegge_FastPi.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/Flatte.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"
#include "AMPTOOLS_AMPS/SinglePS.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF0.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF2.h"
#include "AMPTOOLS_AMPS/KopfKMatrixA0.h"
#include "AMPTOOLS_AMPS/KopfKMatrixA2.h"
#include "AMPTOOLS_AMPS/KopfKMatrixRho.h"
#include "AMPTOOLS_AMPS/KopfKMatrixPi1.h"
#include "AMPTOOLS_AMPS/Vec_ps_moment.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpToolsMPI/AmpToolsInterfaceMPI.h"
#include "IUAmpToolsMPI/DataReaderMPI.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

int rank_mpi;
int size;

double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile) {
   AmpToolsInterfaceMPI ati( cfgInfo );
   bool fitFailed = true;
   double lh = 1e7;

   if(rank_mpi==0) {
      cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

      MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
      fitManager->setMaxIterations(maxIter);

      if( useMinos )
         fitManager->minosMinimization();
      else
         fitManager->migradMinimization();

      if(hesse)
         fitManager->hesseEvaluation();

      fitFailed = ( fitManager->status() != 0 || fitManager->eMatrixStatus() != 3 );

      if( fitFailed )
         cout << "ERROR: fit failed use results with caution..." << endl;
      else 
         lh = ati.likelihood();

      cout << "LIKELIHOOD AFTER MINIMIZATION:  " << lh << endl;
      ati.finalizeFit();
   }


   if( rank_mpi==0 && seedfile.size() != 0 && !fitFailed )
      ati.fitResults()->writeSeed( seedfile );

   ati.exitMPI();
   MPI_Finalize();

   return lh;
}

void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, int numRnd, double maxFraction) {
   AmpToolsInterfaceMPI ati( cfgInfo );

   MinuitMinimizationManager* fitManager = NULL; 
   vector< vector<string> > parRangeKeywords;
   double minLH;
   int minFitTag;
   string fitName;

   if(rank_mpi==0) {
      fitName = cfgInfo->fitName();

      cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

      fitManager = ati.minuitMinimizationManager();
      fitManager->setMaxIterations(maxIter);

      parRangeKeywords = cfgInfo->userKeywordArguments("parRange");

      // keep track of best fit (mininum log-likelihood)
      minLH = numeric_limits<double>::max();
      minFitTag = -1;
   }

   for(int i=0; i<numRnd; i++) {
      bool fitFailed = true;
      double curLH = 1e7;

      if(rank_mpi==0) {
         cout << endl << "###############################" << endl;
         cout << "FIT " << i << " OF " << numRnd << endl;
         cout << endl << "###############################" << endl;

	 // re-initialize parameters from configuration file (reset those not randomized)
	 ati.reinitializePars();

         // randomize parameters
         ati.randomizeProductionPars(maxFraction);
         for(size_t ipar=0; ipar<parRangeKeywords.size(); ipar++) {
            ati.randomizeParameter(parRangeKeywords[ipar][0], atof(parRangeKeywords[ipar][1].c_str()), atof(parRangeKeywords[ipar][2].c_str()));
         }

         if(useMinos)
            fitManager->minosMinimization();
         else
            fitManager->migradMinimization();

         if(hesse)
            fitManager->hesseEvaluation();

         fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

         if( fitFailed )
            cout << "ERROR: fit failed use results with caution..." << endl;

         curLH = ati.likelihood();
         cout << "LIKELIHOOD AFTER MINIMIZATION:  " << curLH << endl;

	 ati.finalizeFit(to_string(i));

         if( seedfile.size() != 0 && !fitFailed ){
            string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
            string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));

            string seedfile_rand = seedfileBaseName + Form("_%d%s", i, seedfileExtension.data());
            ati.fitResults()->writeSeed( seedfile_rand );
         }

         // update best fit
         if( !fitFailed && curLH < minLH ) {
            minLH = curLH;
            minFitTag = i;
         }
      }
   }


   // print best fit results
   if(rank_mpi==0) {
      if(minFitTag < 0) cout << "ALL FITS FAILED!" << endl;
      else {
         cout << "MINIMUM LIKELIHOOD FROM " << minFitTag << " of " << numRnd << " RANDOM PRODUCTION PARS = " << minLH << endl;
         gSystem->Exec(Form("cp %s_%d.fit %s.fit", fitName.data(), minFitTag, fitName.data()));
         if( seedfile.size() != 0 ) {
            string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
            string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));
            gSystem->Exec(Form("cp %s_%d%s %s", seedfileBaseName.data(), minFitTag, seedfileExtension.data(), seedfile.data()));
         }
      }
   }

   ati.exitMPI();
   MPI_Finalize();
}


void runParScan(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, string parScan) {
   double minVal=0, maxVal=0, stepSize=0;
   int steps=0;

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

   AmpToolsInterfaceMPI ati( cfgInfo );
   string fitName;
   ParameterManager* parMgr = NULL;
   MinuitMinimizationManager* fitManager = NULL;

   if(rank_mpi==0) {
      fitName = cfgInfo->fitName();
      cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

      parMgr = ati.parameterManager();
      fitManager = ati.minuitMinimizationManager();
      fitManager->setMaxIterations(maxIter);
   }

   for(int i=0; i<steps; i++) {
      bool fitFailed = true;
      double curLH = 1e7;

      if(rank_mpi==0) {
         cout << endl << "###############################" << endl;
         cout << "FIT " << i << " OF " << steps << endl;
         cout << endl << "###############################" << endl;

	 // reinitialize production parameters from seed file
	 ati.reinitializePars();

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
         cfgInfo->setFitName( fitName + "_scan" );

         if(useMinos)
            fitManager->minosMinimization();
         else
            fitManager->migradMinimization();

         if(hesse)
            fitManager->hesseEvaluation();

         fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);
         curLH = ati.likelihood();

         if( fitFailed )
            cout << "ERROR: fit failed use results with caution..." << endl;

         cout << "LIKELIHOOD AFTER MINIMIZATION:  " << curLH << endl;

         ati.finalizeFit(to_string(i));
         if( seedfile.size() != 0 && !fitFailed ){
            string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
            string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));

            string seedfile_scan = seedfileBaseName + Form("_scan_%d%s", i, seedfileExtension.data());
            ati.fitResults()->writeSeed( seedfile_scan );
         }
      }
   }
   ati.exitMPI();
   MPI_Finalize();
}

void getLikelihood( ConfigurationInfo* cfgInfo ){
    AmpToolsInterfaceMPI ati( cfgInfo );
    if( rank_mpi == 0 )
        cout << "LIKELIHOOD WITHOUT MINIMIZATION:  " << ati.likelihood() << endl;
    ati.exitMPI();
    MPI_Finalize();
}

int main( int argc, char* argv[] ){

   MPI_Init( &argc, &argv );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank_mpi );
   MPI_Comm_size( MPI_COMM_WORLD, &size );

   // set default parameters

   bool useMinos = false;
   bool hesse = false;
   bool noFit = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   unsigned int randomSeed=static_cast<unsigned int>(time(NULL));
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
      if (arg == "-rs"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  randomSeed = atoi(argv[++i]); }
      if (arg == "-m"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  maxIter = atoi(argv[++i]); }
      if (arg == "-n") useMinos = true;
      if (arg == "-H") hesse = true;
      if (arg == "-l") noFit = true;
      if (arg == "-p"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  scanPar = argv[++i]; }
      if (arg == "-h"){
         if(rank_mpi==0) {
            cout << endl << " Usage for: " << argv[0] << endl << endl;
            cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
            cout << "   -H \t\t\t\t\t evaluate HESSE matrix after minimization" << endl;
            cout << "   -c <file>\t\t\t\t config file" << endl;
            cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
            cout << "   -r <int>\t\t\t Perform <int> fits each seeded with random parameters" << endl;
            cout << "   -rs <int>\t\t\t Sets the random seed used by the random number generator for the fits with randomized initial parameters. If not set will use the time()" << endl;
            cout << "   -p <parameter> \t\t\t\t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
            cout << "   -m <int>\t\t\t Maximum number of fit iterations" << endl; 
            cout << "   -l \t\t\t\t Calculate likelihood and exit without running a fit" << endl; 
         }
         MPI_Finalize();
         exit(1);
      }
   }

   if (configfile.size() == 0){
      cout << "No config file specified" << endl;
      MPI_Finalize();
      exit(1);
   }

   ConfigFileParser parser(configfile);
   ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
   if( rank_mpi == 0 ) cfgInfo->display();

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
   AmpToolsInterface::registerAmplitude( VecRadiative_SDME() );
   AmpToolsInterface::registerAmplitude( Zlm() );
   AmpToolsInterface::registerAmplitude( b1piAngAmp() );
   AmpToolsInterface::registerAmplitude( polCoef() );
   AmpToolsInterface::registerAmplitude( Uniform() );
   AmpToolsInterface::registerAmplitude( DblRegge_FastEta() );
   AmpToolsInterface::registerAmplitude( DblRegge_FastPi() );
   AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
   AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
   AmpToolsInterface::registerAmplitude( Piecewise() );
   AmpToolsInterface::registerAmplitude( Flatte() );
   AmpToolsInterface::registerAmplitude( PhaseOffset() );
   AmpToolsInterface::registerAmplitude( ComplexCoeff() );
   AmpToolsInterface::registerAmplitude( OmegaDalitz() );
   AmpToolsInterface::registerAmplitude( LowerVertexDelta() );
   AmpToolsInterface::registerAmplitude( SinglePS() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixF0() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixF2() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixA0() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixA2() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixRho() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixPi1() );
   AmpToolsInterface::registerAmplitude( Vec_ps_moment() );

   AmpToolsInterface::registerDataReader( DataReaderMPI<ROOTDataReader>() );
   AmpToolsInterface::registerDataReader( DataReaderMPI<ROOTDataReaderBootstrap>() );
   AmpToolsInterface::registerDataReader( DataReaderMPI<ROOTDataReaderWithTCut>() );
   AmpToolsInterface::registerDataReader( DataReaderMPI<ROOTDataReaderTEM>() );
   AmpToolsInterface::registerDataReader( DataReaderMPI<FSRootDataReader>() );

   // normalize seedFile file extension if not already set by user
   if (seedfile.size() != 0){
      if (seedfile.find(".") == string::npos) seedfile += ".txt"; 
      else{
         // this ensures the file has an extension for cases like "./my_seed_file" 
         string fileExtension = seedfile.substr(seedfile.find_last_of("."));
         if (fileExtension.find("/") != string::npos) seedfile += ".txt";
      }
   }

   if(noFit)
      getLikelihood(cfgInfo);
   else if(numRnd==0){
      if(scanPar=="")
         runSingleFit(cfgInfo, useMinos, hesse, maxIter, seedfile);
      else
         runParScan(cfgInfo, useMinos, hesse, maxIter, seedfile, scanPar);
   } else {
      cout << "Running " << numRnd << " fits with randomized parameters with seed=" << randomSeed << endl;
      AmpToolsInterface::setRandomSeed(randomSeed);
      runRndFits(cfgInfo, useMinos, hesse, maxIter, seedfile, numRnd, 0.5);
   }

   return 0;
}