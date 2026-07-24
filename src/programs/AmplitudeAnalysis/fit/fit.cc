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
#include "AMPTOOLS_DATAIO/FSRootDataReaderBootstrap.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_DoubleSDMEs.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_factorized.h"
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
#include "AMPTOOLS_AMPS/PiPiSWaveAMPK.h"
#include "AMPTOOLS_AMPS/b1piAngAmp.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/polCoef.h"
#include "AMPTOOLS_AMPS/DblRegge_FastEta.h"
#include "AMPTOOLS_AMPS/DblRegge_FastPi.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/Iso_ps_refl.h"
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

void summarizeFits(std::vector<std::tuple<int,bool,int,int,double>>& fitLLs, bool csvOutput = false,
                    double tieEps = 1e-3) {
    // tieEps: LogL difference below which two converged fits are considered "tied".
    // MINUIT's default EDM convergence criterion is on the order of 1e-3 (for UP=1,
    // i.e. -logL directly), so two independent runs that land on the same true
    // minimum will typically agree to within about this tolerance. Tighten this
    // (e.g. 1e-4) if you called MIGRAD with a smaller "tolerance" arg / SetErrorDef,
    // or loosen it if you want a "statistically indistinguishable" cut instead
    // (e.g. ~2.0, akin to a likelihood-ratio / AIC-difference threshold).

    if(fitLLs.empty()) return;

    // sort: converged (success) fits first, ranked by LogL ascending;
    // failed fits pushed to the bottom (also LogL-sorted among themselves, but unranked)
    std::sort(fitLLs.begin(), fitLLs.end(),
        [](const auto& a, const auto& b){
            bool sa = std::get<1>(a), sb = std::get<1>(b);
            if (sa != sb) return sa; // successes before failures
            return std::get<4>(a) < std::get<4>(b);
        });

    // collect all *successful* fits sharing the minimum LogL (ties, within tieEps)
    auto tie_min = [&](int idx) {
        using FitType = std::tuple<int,bool,int,int,double>;
        double minVal = std::numeric_limits<double>::max();
        for (const auto& v : fitLLs) {
            if (!std::get<1>(v)) continue; // only consider converged fits
            double val = (idx == 4) ? std::get<4>(v) : (double)std::get<3>(v);
            if (val < minVal) minVal = val;
        }
        std::vector<FitType> out;
        for (const auto& v : fitLLs) {
            if (!std::get<1>(v)) continue;
            double val = (idx == 4) ? std::get<4>(v) : (double)std::get<3>(v);
            if (std::fabs(val - minVal) < tieEps) out.push_back(v);
        }
        return out;
    };

    // success rate
    size_t succ = 0;
    for (const auto& f : fitLLs) if (std::get<1>(f)) ++succ;
    double pct = 100.0 * static_cast<double>(succ) / static_cast<double>(fitLLs.size());

    // best successful fit (lowest LogL among successes only)
    int    bestSuccIdx  = -1;
    double bestSuccLogL = std::numeric_limits<double>::max();
    for (const auto& f : fitLLs) {
        if (std::get<1>(f) && std::get<4>(f) < bestSuccLogL) {
            bestSuccLogL = std::get<4>(f);
            bestSuccIdx  = std::get<0>(f);
        }
    }

    auto bestByLogL = tie_min(4);

    // reference LogL for DeltaLogL column: best converged fit if one exists,
    // otherwise just fall back to the first row so the column still prints
    double refLogL = (bestSuccIdx >= 0) ? bestSuccLogL : std::get<4>(fitLLs[0]);

    // ---- write to both stdout and file ----
    std::string outfile = csvOutput ? "fit_ranking.csv" : "fit_ranking.txt";
    std::ofstream fout(outfile);
    if (!fout) std::cerr << "Error: cannot open " << outfile << "\n";

    auto write = [&](std::ostream& os) {

        if (csvOutput) {
            os << "Rank,FitIndex,Success,FitStatus,eMatrix,LogL,DeltaLogL\n";
            int rank = 0;
            for (size_t i = 0; i < fitLLs.size(); ++i) {
                const auto& f = fitLLs[i];
                bool isSucc = std::get<1>(f);
                if (isSucc) ++rank;
                os << (isSucc ? std::to_string(rank) : std::string("-")) << ","
                   << std::get<0>(f)          << ","
                   << (isSucc ? 1:0)          << ","
                   << std::get<2>(f)          << ","
                   << std::get<3>(f)          << ","
                   << std::setprecision(10) << std::fixed << std::get<4>(f) << ","
                   << std::setprecision(6)  << std::fixed << (std::get<4>(f) - refLogL)
                   << "\n";
            }
        } else {
            os << "\n===== SUMMARY OF ALL FITS (converged fits ranked by LogL; failed fits listed unranked) =====\n";
            os << std::left
               << std::setw(6)  << "Rank"
               << std::setw(6)  << "#"
               << std::setw(9)  << "Success"
               << std::setw(11) << "FitStatus"
               << std::setw(9)  << "eMatrix"
               << std::setw(18) << "LogL"
               << std::setw(14) << "DeltaLogL"
               << "\n";
            int rank = 0;
            for (size_t i = 0; i < fitLLs.size(); ++i) {
                const auto& f = fitLLs[i];
                bool isSucc = std::get<1>(f);
                if (isSucc) ++rank;
                std::string rankStr = isSucc ? std::to_string(rank) : std::string("-");
                os << std::left
                   << std::setw(6)  << rankStr
                   << std::setw(6)  << std::get<0>(f)
                   << std::setw(9)  << (isSucc ? "Y" : "N")
                   << std::setw(11) << std::get<2>(f)
                   << std::setw(9)  << std::get<3>(f)
                   << std::setw(18) << (std::ostringstream() << std::setprecision(10) << std::fixed << std::get<4>(f)).str()
                   << std::setw(14) << (std::ostringstream() << std::setprecision(6)  << std::fixed << (std::get<4>(f) - refLogL)).str()
                   << "\n";
            }
        }

        // best fit callout (converged fits only)
        if (!bestByLogL.empty()) {
            os << "\nBEST FIT(S) by LogL among converged fits (" << bestByLogL.size() << " tied):\n";
            for (const auto& f : bestByLogL)
                os << "  Fit #" << std::get<0>(f)
                   << "  success=" << (std::get<1>(f) ? "Y":"N")
                   << "  status="  << std::get<2>(f)
                   << "  eMatrix=" << std::get<3>(f)
                   << "  LogL="    << std::setprecision(10) << std::fixed << std::get<4>(f)
                   << "\n";
        } else {
            os << "\nBEST FIT(S) by LogL among converged fits: none (no fits converged)\n";
        }

        if (bestSuccIdx >= 0)
            os << "BEST SUCCESSFUL FIT: #" << bestSuccIdx
               << "  LogL=" << std::setprecision(10) << std::fixed << bestSuccLogL << "\n";
        else
            os << "BEST SUCCESSFUL FIT: none (all fits failed)\n";

        // success rate
        os << "\nSuccess: " << succ << " / " << fitLLs.size()
           << " (" << std::setprecision(1) << std::fixed << pct << "%)\n";
        os << "=================================================\n";
    };

    write(std::cout);
    if (fout) write(fout);
}

double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, int eMatrixRequirement) {
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
    ( fitManager->status() != 0 || fitManager->eMatrixStatus() < eMatrixRequirement );

  if( fitFailed ){
    cout << "ERROR: fit failed use results with caution..." << endl;
    cout << "Fit status = " << fitManager->status() << ", eMatrix status = " << fitManager->eMatrixStatus() << endl;
    return 1e6;
  }

  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

  ati.finalizeFit();

  if( seedfile.size() != 0 && !fitFailed ){
    ati.fitResults()->writeSeed( seedfile );
  }

  return ati.likelihood();
}

void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, int numRnd, double maxFraction, int eMatrixRequirement) {
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

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() < eMatrixRequirement);

    if( fitFailed ){
      cout << "ERROR: fit failed use results with caution..." << endl;
      cout << "Fit status = " << fitManager->status() << ", eMatrix status = " << fitManager->eMatrixStatus() << endl;
    }

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
    // summarizeFits(fitLLs,true); // write summary to CSV file
    summarizeFits(fitLLs);
  }
}

void runParScan(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, string parScan, int eMatrixRequirement) {
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

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() < eMatrixRequirement);

    if( fitFailed ){
      cout << "ERROR: fit failed use results with caution..." << endl;
      cout << "Fit status = " << fitManager->status() << ", eMatrix status = " << fitManager->eMatrixStatus() << endl;
    }

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

void saveDummyFit( ConfigurationInfo* cfgInfo ){
  AmpToolsInterface ati( cfgInfo );
  ati.finalizeFit();
  return;
}


int main( int argc, char* argv[] ){

   // set default parameters

   bool useMinos = false;
   bool hesse = false;
   bool noFit = false;
   bool printAmps = false;
   bool createDummyFit = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   unsigned int randomSeed=static_cast<unsigned int>(time(NULL));
   int maxIter = 10000;
   int eMatrixRequirement = 3;

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
      if (arg == "-e"){
         if ((i+1 ==argc)  || (argv[i+1][0] == '-')) arg = "-h";
         else  eMatrixRequirement = atoi(argv[++i]); }
      if (arg == "-n") useMinos = true;
      if (arg == "-H") hesse = true;
      if (arg == "-l") noFit = true;
      if (arg == "-test" ) printAmps = true;
      if (arg == "-d") createDummyFit = true;
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
         cout << "   -p <parameter> \t\t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
         cout << "   -m <int>\t\t\t Maximum number of fit iterations" << endl; 
         cout << "   -e <int>\t\t\t Minimum required level of error matrix status for a successful fit." << endl;
         cout << "   \t\t\t\t\t 0 = not calculated at all" << endl;
         cout << "   \t\t\t\t\t 1 = approximation only, not accurate" << endl;
         cout << "   \t\t\t\t\t 2 = full matrix, but forced positive-definite" << endl;
         cout << "   \t\t\t\t\t 3 = full accurate covariance matrix (default, recommended for most fits)" << endl;
         cout << "   -l \t\t\t\t Calculate likelihood and exit without running a fit" << endl; 
	       cout << "   -test \t\t\t Print amplitude details for the first 2 data events" << endl;
         cout << "   -d \t\t\t\t Create dummy .fit file and exit without running a fit." << endl;         
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
   AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_DoubleSDMEs());
   AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_factorized());
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
   AmpToolsInterface::registerAmplitude( Iso_ps_refl() );
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
   AmpToolsInterface::registerAmplitude( PiPiSWaveAMPK() );
   
   AmpToolsInterface::registerDataReader( ROOTDataReader() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderHist() );
   AmpToolsInterface::registerDataReader( FSRootDataReader() );
   AmpToolsInterface::registerDataReader( FSRootDataReaderBootstrap() );

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
   else if(createDummyFit)
     saveDummyFit(cfgInfo);
   else if(printAmps)
     printAmplitudes(cfgInfo);
   else if(numRnd==0){
     if(scanPar=="")
       runSingleFit(cfgInfo, useMinos, hesse, maxIter, seedfile, eMatrixRequirement);
     else
       runParScan(cfgInfo, useMinos, hesse, maxIter, seedfile, scanPar, eMatrixRequirement);
   } else {
     cout << "Running " << numRnd << " fits with randomized parameters with seed=" << randomSeed << endl;
     AmpToolsInterface::setRandomSeed(randomSeed);
     runRndFits(cfgInfo, useMinos, hesse, maxIter, seedfile, numRnd, 0.5, eMatrixRequirement);
   }

  return 0;
}


