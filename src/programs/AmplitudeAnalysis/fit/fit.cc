#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <limits>
#include <tuple>   

#include "TSystem.h"

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderWithTCut.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderHist.h"
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
#include "AMPTOOLS_AMPS/TwoLeptonAngles.h"
#include "AMPTOOLS_AMPS/TwoLeptonAnglesGJ.h"
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
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"
#include "AMPTOOLS_AMPS/SinglePS.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF0.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF2.h"
#include "AMPTOOLS_AMPS/KopfKMatrixA0.h"
#include "AMPTOOLS_AMPS/KopfKMatrixA2.h"
#include "AMPTOOLS_AMPS/KopfKMatrixRho.h"
#include "AMPTOOLS_AMPS/KopfKMatrixPi1.h"
#include "AMPTOOLS_AMPS/Vec_ps_moment.h"
#include "UTILITIES/randomized_sdme.h"


#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

void summarizeFits(std::vector<std::tuple<int,bool,int,int,double>>& fitLLs) {
    if(fitLLs.size() == 0) return;
    std::sort(fitLLs.begin(), fitLLs.end(), [](const std::tuple<int,bool,int,int,double>& a, const std::tuple<int,bool,int,int,double>& b){ return std::get<4>(a) < std::get<4>(b); });
    std::ofstream fout("fit_ranking.txt");
    if(!fout) std::cerr << "Error: cannot open fit_ranking.txt\n";
    auto print_header = [](std::ostream& os){ os << "\nSUMMARY OF ALL FITS:\n" << std::left << std::setw(6) << "#" << std::setw(9) << "Success" << std::setw(11) << "FitStatus" << std::setw(9) << "eMatrix" << std::setw(14) << "LogL" << '\n'; };
    auto print_row = [&](std::ostream& os, size_t i){ os << std::left << std::setw(6) << std::get<0>(fitLLs[i]) << std::setw(9) << (std::get<1>(fitLLs[i]) ? "Y" : "N") << std::setw(11) << std::get<2>(fitLLs[i]) << std::setw(9) << std::get<3>(fitLLs[i]) << std::setw(14) << (std::ostringstream() << std::setprecision(8) << std::defaultfloat << std::get<4>(fitLLs[i])).str() << '\n'; };
    print_header(std::cout);
    if(fout) print_header(fout);
    for(size_t i=0; i<fitLLs.size(); ++i) {
        print_row(std::cout, i);
        if(fout) print_row(fout, i);
    }
    size_t succ = 0; // print overall success rate
    for (size_t i = 0; i < fitLLs.size(); ++i) if (std::get<1>(fitLLs[i])) ++succ;
    double pct = fitLLs.empty() ? 0.0 : (100.0 * static_cast<double>(succ) / static_cast<double>(fitLLs.size()));
    std::cout << "\nSuccess: " << succ << " / " << fitLLs.size() << " (" << std::setprecision(3) << std::fixed << pct << "%)\n";
}

double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile) {
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

  if( hesse )
     fitManager->hesseEvaluation();

  bool fitFailed =
    ( fitManager->status() != 0 || fitManager->eMatrixStatus() != 3 );

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

void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, int numRnd, double maxFraction) {
  AmpToolsInterface ati( cfgInfo );
  string fitName = cfgInfo->fitName();

  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);

  vector< vector<string> > parRangeKeywords = cfgInfo->userKeywordArguments("parRange");
  vector< vector<string> > parSDMEKeywords = cfgInfo->userKeywordArguments("parSDME");


  // keep track of best fit (mininum log-likelihood)
  double minLL = numeric_limits<double>::max();
  int minFitTag = -1;

  // tuple <index, success_flag, fit_status, eMatrix_status, loglikelihood>
  vector < tuple<int, bool, int, int, double> > fitLLs; 

  for(int i=0; i<numRnd; i++) {
    cout << endl << "###############################" << endl;
    cout << "FIT " << i << " OF " << numRnd << endl;
    cout << endl << "###############################" << endl;

    // re-initialize parameters from configuration file (reset those not randomized)
    ati.reinitializePars();

    // set maximal fraction of production parameters to randomize
    cout << "Maximal fraction allowed: " << std::setprecision(2) << std::defaultfloat << maxFraction << endl;
    ati.randomizeProductionPars(maxFraction);

    // Initiate SDME values from randomized hermitian matrices  
    if (parSDMEKeywords.size() > 0) {
      assert(parSDMEKeywords.size() == 9 && "randomized_sdme() returns exactly 9 SDMEs");
      std::cout << "Initializing " << parSDMEKeywords.size() << " SDME values from randomized hermitian matrices..." << std::endl;
      // Order: rho000, rho100, rho1m10, rho111, rho001, rho101, rho1m11, rho102, rho1m12
      std::vector<double> sdme_values = randomized_sdme(true); // use normal distribution by default, change to false for uniform distribution
      for (size_t ipar = 0; ipar < parSDMEKeywords.size(); ipar++) {
        double init_value = sdme_values[ipar];
        std::cout << "Initializing " << std::setprecision(3) << std::defaultfloat << parSDMEKeywords[ipar][0] << " to value " << init_value << std::endl;
        ati.randomizeParameter(parSDMEKeywords[ipar][0], init_value, init_value);
      }
    }

    // randomize other user-defined parameters
    if (parRangeKeywords.size() > 0){ 
      std::cout << "Initializing parRange parameters... " << std::endl;
      std::cout << "Number of user parameters to randomize: " << parRangeKeywords.size() << std::endl;
      for (size_t ipar = 0; ipar < parRangeKeywords.size(); ipar++) {
        std::cout << "Randomizing parameter: " 
        << parRangeKeywords[ipar][0] << " uniformly in range [" 
        << parRangeKeywords[ipar][1] << ", " 
        << parRangeKeywords[ipar][2] << "]" << std::endl;
        ati.randomizeParameter(parRangeKeywords[ipar][0], atof(parRangeKeywords[ipar][1].c_str()), atof(parRangeKeywords[ipar][2].c_str()));
      }
    }

    if(useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if(hesse)
       fitManager->hesseEvaluation();

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

    if( fitFailed )
      cout << "ERROR: fit failed use results with caution..." << endl;

    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

    ati.finalizeFit(to_string(i));
    fitLLs.push_back(std::make_tuple(
        i, !fitFailed, fitManager->status(), fitManager->eMatrixStatus(), ati.likelihood()
    ));

    if( seedfile.size() != 0 && !fitFailed ){
      string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
      string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));

      string seedfile_rand = seedfileBaseName + Form("_%d%s", i, seedfileExtension.data());
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
    if( seedfile.size() != 0 ){
      string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
      string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));
      gSystem->Exec(Form("cp %s_%d%s %s", seedfileBaseName.data(), minFitTag, seedfileExtension.data(), seedfile.data()));
    }
    // print summary of all fits
    summarizeFits(fitLLs);
  }
}

void runBootstrapFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, int numBootstrap, unsigned int bootstrapSeed) {
  /* 
  I can actually check what dataReaders the cfgInfo data files are read with, and 
  replace it with the bootstrap one if not set. 
  
  Trouble is that the AmpToolsInterface will always run resetConfigurationInfo if
  we pass a cfgInfo object to it. Somewhere in here, this will run the dataReaders and 
  thus load in the data again. The loadEvents explicitly will call the data reader
  loading. If I can clear and reload just the data events with a new dataReader, this 
  would be ideal. I'll need to debug to track down when this happens
  
  I would effectively be doing the same thing if I create a new AmpToolsInterface object
  for each bootstrap fit. Somehow I need to load the interface without it trying to load
  in the new data.

  Also, I should consider allowing a seedfile to be passed here. Would make it easier
  for user to simply give this the same fit.cfg and seedfile from the rand results, and
  for this function to handle that for them. All I have to do is have the seedfile be
  included at the cfgInfo parser, which might have an option to access that directly
  
  */
  AmpToolsInterface ati( cfgInfo );

  std::cout << "first likelihood "<< ati.likelihood() << std::endl;

  return ati.likelihood();

  // ReactionInfo* reaction = cfgInfo->reaction("omegapi");
  // std::string oldDataFile;
  // std::string oldDataReader;
  // oldDataReader = reaction->data().first;
  // oldDataFile   = reaction->data().second[0];

  // std::vector<std::string> new_data = {"anglesOmegaPiAmplitude_45.root"};
  // reaction->setData(oldDataReader, new_data);  

  // std::cout << "Replaced data, call likelihood before reset: " << ati.likelihood() << std::endl;

  // cout << reaction->genMC().first << " " << reaction->genMC().second[0] << endl;
  // cout << reaction->accMC().first << " " << reaction->accMC().second[0] << endl;

  // ReactionInfo* new_reaction = cfgInfo->reaction("omegapi");
  // cout << new_reaction->genMC().first << " " << new_reaction->genMC().second[0] << endl;
  // cout << new_reaction->accMC().first << " " << new_reaction->accMC().second[0] << endl;

  // ati.resetConfigurationInfo(cfgInfo);

  // std::cout << "Replaced data, call likelihood after reset: " << ati.likelihood() << std::endl;
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

    cfgInfo->setFitName(fitName + "_scan");

    if(useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if(hesse)
       fitManager->hesseEvaluation();

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

    if( fitFailed )
      cout << "ERROR: fit failed use results with caution..." << endl;

    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

    ati.finalizeFit(to_string(i));

    if( seedfile.size() != 0 && !fitFailed ){
      string seedfileBaseName = seedfile.substr(0, seedfile.find_last_of("."));
      string seedfileExtension = seedfile.substr(seedfile.find_last_of("."));

      string seedfile_scan = seedfileBaseName + Form("_scan_%d%s", i, seedfileExtension.data());
      ati.fitResults()->writeSeed( seedfile_scan );
    }
  }
}

void getLikelihood( ConfigurationInfo* cfgInfo ){
    AmpToolsInterface ati( cfgInfo );
    cout << "LIKELIHOOD WITHOUT MINIMIZATION:  " << ati.likelihood() << endl;
    return;
}

void printAmplitudes( ConfigurationInfo* cfgInfo ){

  AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kPlotGeneration );
  DataReader* dataReader = ati.dataReader(cfgInfo->reactionList()[0]->reactionName());  
  dataReader->resetSource();
  for (int i = 0; i < 2; i++){
    Kinematics* kin = dataReader->getEvent();
    ati.printEventDetails(cfgInfo->reactionList()[0]->reactionName(),kin);
    delete kin;
  }
  return;
}


int main( int argc, char* argv[] ){

   // set default parameters

   bool useMinos = false;
   bool hesse = false;
   bool noFit = false;
   bool printAmps = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   int numBootstrap = 0;
   unsigned int randomSeed=static_cast<unsigned int>(time(NULL));
   unsigned int bootstrapSeed=static_cast<unsigned int>(time(NULL)); // TODO: implement this, for now not used since bootstrapreader requires an int
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
      if (arg == "-b"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  numBootstrap = atoi(argv[++i]); }
      if (arg == "-bs"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  bootstrapSeed = atoi(argv[++i]); }
      if (arg == "-m"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  maxIter = atoi(argv[++i]); }
      if (arg == "-n") useMinos = true;
      if (arg == "-H") hesse = true;
      if (arg == "-l") noFit = true;
      if (arg == "-test" ) printAmps = true;
      if (arg == "-p"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  scanPar = argv[++i]; }
      if (arg == "-h"){
         cout << endl << " Usage for: " << argv[0] << endl << endl;
         cout << "   -n \t\t\t\t use MINOS instead of MIGRAD" << endl;
         cout << "   -H \t\t\t\t evaluate HESSE matrix after minimization" << endl;
         cout << "   -c <file>\t\t\t config file" << endl;
         cout << "   -s <output file>\t\t for seeding next fit based on this fit (optional)" << endl;
         cout << "   -r <int>\t\t\t Perform <int> fits each seeded with random parameters" << endl;
         cout << "   -rs <int>\t\t\t Sets the random seed used by the random number generator for the fits with randomized initial parameters. If not set will use the time()" << endl;
         cout << "   -b <int>\t\t\t Perform <int> bootstrap fits" << endl;
         cout << "   -bs <int>\t\t\t Sets the first random seed used by the random number generator for the bootstrap fits. Each subsequent bootstrap fit will use this seed + the bootstrap iteration number. If not set will use the time()" << endl;
         cout << "   -p <parameter> \t\t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
         cout << "   -m <int>\t\t\t Maximum number of fit iterations" << endl; 
         cout << "   -l \t\t\t\t Calculate likelihood and exit without running a fit" << endl; 
	 cout << "   -test \t\t\t Print amplitude details for the first 2 data events" << endl;
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
   AmpToolsInterface::registerAmplitude( TwoPiNC_tdist() );
   AmpToolsInterface::registerAmplitude( ThreePiAngles() );
   AmpToolsInterface::registerAmplitude( ThreePiAnglesSchilling() );
   AmpToolsInterface::registerAmplitude( VecRadiative_SDME() );
   AmpToolsInterface::registerAmplitude( TwoLeptonAngles() );
   AmpToolsInterface::registerAmplitude( TwoLeptonAnglesGJ() );
   AmpToolsInterface::registerAmplitude( Zlm() );
   AmpToolsInterface::registerAmplitude( b1piAngAmp() );
   AmpToolsInterface::registerAmplitude( polCoef() );
   AmpToolsInterface::registerAmplitude( Uniform() );
   AmpToolsInterface::registerAmplitude( DblRegge_FastEta() );
   AmpToolsInterface::registerAmplitude( DblRegge_FastPi() );
   AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
   AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
   AmpToolsInterface::registerAmplitude( PhaseOffset() );
   AmpToolsInterface::registerAmplitude( ComplexCoeff() );
   AmpToolsInterface::registerAmplitude( OmegaDalitz() );
   AmpToolsInterface::registerAmplitude( Piecewise() );
   AmpToolsInterface::registerAmplitude( LowerVertexDelta() );
   AmpToolsInterface::registerAmplitude( SinglePS() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixF0() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixF2() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixA0() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixA2() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixRho() );
   AmpToolsInterface::registerAmplitude( KopfKMatrixPi1() );
   AmpToolsInterface::registerAmplitude( Vec_ps_moment() );

   AmpToolsInterface::registerDataReader( ROOTDataReader() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderHist() );
   AmpToolsInterface::registerDataReader( FSRootDataReader() );

   // normalize seedFile file extension if not already set by user
   if (seedfile.size() != 0){
     if (seedfile.find(".") == string::npos) seedfile += ".txt";
     else{
       // this ensures the file has an extension for cases like "./my_seed_file" 
       string fileExtension = seedfile.substr(seedfile.find_last_of("."));
       if (fileExtension.find("/") != string::npos) seedfile += ".txt";
     }
   }

   if(numRnd > 0 && numBootstrap > 0) {
     cout << "Error: cannot perform both random parameter fits and bootstrap fits simultaneously. Please pick one or the other. Aborting." << endl;
     return 1;
   }

   if(noFit)
     getLikelihood(cfgInfo);
   else if(printAmps)
     printAmplitudes(cfgInfo);
   else if(numRnd==0 && numBootstrap==0){
     if(scanPar=="")
       runSingleFit(cfgInfo, useMinos, hesse, maxIter, seedfile);
     else
       runParScan(cfgInfo, useMinos, hesse, maxIter, seedfile, scanPar);
   } else if(numBootstrap > 0){
       cout << "Running " << numBootstrap << " bootstrap fits beginning at seed=" << bootstrapSeed << endl;
       runBootstrapFits(cfgInfo, useMinos, hesse, maxIter, numBootstrap, bootstrapSeed);
   } else {
     cout << "Running " << numRnd << " fits with randomized parameters with seed=" << randomSeed << endl;
     AmpToolsInterface::setRandomSeed(randomSeed);
     runRndFits(cfgInfo, useMinos, hesse, maxIter, seedfile, numRnd, 0.5);
   }

  return 0;
}


