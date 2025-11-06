#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <limits>
#include <cmath>
#include <random>
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
#include "AMPTOOLS_AMPS/TwoPSMoment.h"
#include "AMPTOOLS_AMPS/VecPSMoment.h"
#include "AMPTOOLS_AMPS/KStarHyperon.h"
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

// std::vector<double> randomized_sdme(bool use_normal_dist = true) {
//     using cplx = std::complex<double>;

//     // Random number generator
//     static std::random_device rd;
//     static std::mt19937 gen(rd());

//     // Distributions
//     static std::normal_distribution<> normal_dist(0.0, 1.0);
//     static std::uniform_real_distribution<> uniform_dist(-1.0, 1.0);

//     // Unified sampling function
//     auto sample = [&]() {
//         return use_normal_dist ? normal_dist(gen) : uniform_dist(gen);
//     };

//     // Step 1: generate 3x3 complex matrix of helicity amplitudes
//     cplx M[3][3];
//     for (int i = 0; i < 3; ++i)
//         for (int j = 0; j < 3; ++j)
//             M[i][j] = cplx(sample(), sample());

//     // Step 2: compute normalization N = 0.5 * sum |M|^2
//     double N = 0.0;
//     for (int i = 0; i < 3; ++i)
//         for (int j = 0; j < 3; ++j)
//             N += std::norm(M[i][j]);
//     N *= 0.5;

//     // Step 3: define shortcuts 
//     const cplx& M_m1_m1 = M[0][0];
//     //const cplx& M_0_m1  = M[1][0];
//     const cplx& M_p1_m1 = M[2][0];

//     const cplx& M_m1_0  = M[0][1];
//     //const cplx& M_0_0   = M[1][1];
//     const cplx& M_p1_0  = M[2][1];

//     const cplx& M_m1_p1 = M[0][2];
//     //const cplx& M_0_p1  = M[1][2];
//     const cplx& M_p1_p1 = M[2][2];

//     // Step 4: compute SDMEs

//     // (B1a)
//     double rho000 = (std::norm(M_p1_0) + std::norm(M_p1_0)) / N;

//     // (B1b)
//     cplx term_b = (M_p1_p1 - M_p1_m1) * std::conj(M_p1_0);
//     double rho100 = 0.5 * term_b.real() / N;

//     // (B1c)
//     double rho1m10 = (M_p1_p1 * std::conj(M_p1_m1)).real() / N;

//     // (B1d)
//     double rho111 = (M_m1_p1 * std::conj(M_p1_p1)).real() / N;

//     // (B1e)
//     double rho001 = (M_m1_0 * std::conj(M_p1_0)).real() / N;

//     // (B1f, B1g)
//     cplx prod_m1p1 = M_m1_p1 * std::conj(M_p1_m1);
//     cplx prod_11   = M_p1_p1 * std::conj(M_m1_m1);
//     double rho1m11 = (prod_m1p1 + prod_11).real() / N;
//     double rho1m12 = (prod_m1p1 - prod_11).imag() / N;

//     // (B1h, B1i)
//     cplx prod_10   = M_m1_p1 * std::conj(M_p1_0);
//     cplx prod_m10  = M_p1_p1 * std::conj(M_m1_0);
//     double rho101 = (prod_10 + prod_m10).real() / N;
//     double rho102 = (prod_10 - prod_m10).imag() / N;

//     // Output vector of 9 SDMEs
//     return {
//         rho000, rho100, rho1m10,
//         rho111, rho001, rho101,
//         rho1m11, rho102, rho1m12
//     };
// }

// std::vector<double> randomized_sdme_new(bool use_normal_dist = true) {
//   // Random number generator
//   static std::random_device rd;
//   static std::mt19937 gen(rd());

//   // Distributions
//   static std::normal_distribution<> normal_dist(0.0, 1.0);
//   static std::uniform_real_distribution<> uniform_dist(-1.0, 1.0);

//   // Unified sampling function
//   auto sample = [&]() {
//       return use_normal_dist ? normal_dist(gen) : uniform_dist(gen);
//   };

//   // generate 2x3 complex matrix:
//   //    helicity amplitudes M_{lambda_g, lambda_V},
//   //    where lambda_g = 1, -1 and lambda_V = 1, 0, -1,
//   //    then compute normalization N = 0.5 * sum |M|^2.
//   std::complex<double> M[2][3];
//   double N = 0.0; 
//   for (int i = 0; i < 2; ++i){
//     for (int j = 0; j < 3; ++j){
//       M[i][j] = std::complex<double>(sample(), sample());
//       N += 0.5 * std::norm(M[i][j]);
//     }
//   }
      

//   // human readable shortcuts 
//   const std::complex<double>& M_m1_m1 = M[0][0];  // M_{-1, -1}
//   const std::complex<double>& M_p1_m1 = M[1][0];  // M_{+1, -1}
//   const std::complex<double>& M_m1_0  = M[0][1];  // M_{-1, 0}
//   const std::complex<double>& M_p1_0  = M[1][1];  // M_{+1, 0}
//   const std::complex<double>& M_m1_p1 = M[0][2];  // M_{-1, +1}
//   const std::complex<double>& M_p1_p1 = M[1][2];  // M_{+1, +1}

//   // compute SDMEs, from [V. MATHIEU et al. PHYS. REV. D 97, 094003 (2018)]
//   double Ninv = 1.0 / (1e-12 + N); // small offset to avoid division by zero
  
//   // (B1a-c)
//   double rho000  = Ninv * (M_p1_0 * std::conj(M_p1_0)).real();               
//   double rho100  = Ninv * ((M_p1_p1 - M_p1_m1) * std::conj(M_p1_0)).real() * 0.5;   
//   double rho1m10 = Ninv * (M_p1_p1 * std::conj(M_p1_m1)).real();          
  
//   // (B1d-e)
//   double rho111  = Ninv * (M_m1_p1 * std::conj(M_p1_p1)).real();
//   double rho001  = Ninv * (M_m1_0 * std::conj(M_p1_0)).real();

//   // (B1fg)
//   std::complex<double> prod_m1p1 = M_m1_p1 * std::conj(M_p1_m1);
//   std::complex<double> prod_11   = M_p1_p1 * std::conj(M_m1_m1);
//   double rho1m11 = Ninv * (prod_m1p1 + prod_11).real();
//   double rho1m12 = Ninv * (prod_m1p1 - prod_11).imag();

//   // (B1hi)
//   double re_prod_10   = (M_m1_p1 * std::conj(M_p1_0)).real();
//   double re_prod_m10  = (M_p1_p1 * std::conj(M_m1_0)).real();
//   double rho101 = Ninv * (re_prod_10 + re_prod_m10);
//   double rho102 = Ninv * (re_prod_10 - re_prod_m10);

//   // Output vector of 9 SDMEs
//   return {
//       rho000, rho100, rho1m10,
//       rho111, rho001, rho101,
//       rho1m11, rho102, rho1m12
//   };
// }



double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, int strategy, string seedfile) {
  AmpToolsInterface ati( cfgInfo );

  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);
  if(strategy != 1) fitManager->setStrategy(strategy);

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


void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, int strategy, string seedfile, int numRnd, double maxFraction, double T) {
  AmpToolsInterface ati(cfgInfo);
  string fitName = cfgInfo->fitName();

  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);
  if (strategy != 1) fitManager->setStrategy(strategy);

  vector< vector<string> > parRangeKeywords = cfgInfo->userKeywordArguments("parRange");
  vector< vector<string> > parSDMEKeywords = cfgInfo->userKeywordArguments("parSDME");



  double minLL = numeric_limits<double>::max();
  int minFitTag = -1;
  vector< pair<int, double> > tagLLs;  // (tagNumber, LL)

  // tuple <index, success_flag, fit_status, eMatrix_status, loglikelihood>
  vector < tuple<int, bool, int, int, double> > fitLLs; 

  for (int i = 0; i < numRnd; i++) {
    cout << endl << "###############################" << endl;
    cout << "FIT " << i << " OF " << numRnd << endl;
    cout << "###############################" << endl;


    ati.reinitializePars();
    cout << "Maximal fraction allowed: " << std::setprecision(2) << std::defaultfloat << maxFraction << endl;
    ati.randomizeProductionPars(maxFraction);

    std::vector<std::string> keywords = cfgInfo->userKeywords();
    std::cout << "User keywords (" << keywords.size() << "):" << std::endl;
    for (const auto& kw : keywords) {
        std::cout << "  - " << kw << std::endl;
    }

    
      
    // Initiate SDME values from randomized hermitian matrices  
    if (parSDMEKeywords.size() > 0) {
      std::cout << "Initializing SDME values from randomized hermitian matrices... " << std::endl;
      std::cout << "Number of parameters to initiate: " << parSDMEKeywords.size() << std::endl;
      std::vector<double> sdme_values = randomized_sdme(true); // use normal distribution by default, change to false for uniform distribution
      for (size_t ipar = 0; ipar < parSDMEKeywords.size(); ipar++) {
        double init_value = sdme_values[ipar];
        std::cout << "Initializing SDME parameter " << parSDMEKeywords[ipar][0] << " to value " << init_value << std::endl;
        ati.randomizeParameter(parSDMEKeywords[ipar][0], init_value, init_value);
      }
    }

    // std::cout << "Number of user parameters to randomize: " << parRangeKeywords.size() << std::endl;
    // std::vector<double> sdme_values = randomized_sdme(false); // false for uniform distribution
    // for (size_t ipar = 0; ipar < parRangeKeywords.size(); ipar++) {
    //     double init_value = sdme_values[ipar];
    //     std::cout << "Initializing parameter " << parRangeKeywords[ipar][0] 
    //               << " to value " << init_value << std::endl;
    //     ati.randomizeParameter(parRangeKeywords[ipar][0], init_value, init_value);
    // }
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
    


    if (useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if (hesse)
      fitManager->hesseEvaluation();

    bool fitFailed = (strategy == 0) ? (fitManager->status() != 0)
                                     : (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

    if (fitFailed)
      cout << "ERROR: fit failed use results with caution..." << endl;

    double LL = ati.likelihood();
    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << LL << endl;

    ati.finalizeFit(to_string(i));
    fitLLs.push_back(std::make_tuple(
      i, !fitFailed, fitManager->status(), fitManager->eMatrixStatus(), ati.likelihood()
    ));

    if (seedfile.size() != 0 && !fitFailed) {
      string seedfile_rand = seedfile + Form("_%d.txt", i);
      ati.fitResults()->writeSeed(seedfile_rand);
    }

    if (!fitFailed) {
      tagLLs.emplace_back(i, LL);
      if (LL < minLL) {
        minLL = LL;
        minFitTag = i;
      }
    }
  }

  if (minFitTag < 0) {
    cout << "ALL FITS FAILED!" << endl;
    return;
  }

  // Remove duplicates by LL (same LL might converge to same point)
  sort(tagLLs.begin(), tagLLs.end(),
  [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
      return a.second < b.second;
  });
  
  vector< pair<int, double> > uniqueLLs;
  for (const auto& p : tagLLs) {
    if (uniqueLLs.empty() || fabs(p.second - uniqueLLs.back().second) > 1e-2) {
      uniqueLLs.push_back(p);
    }
  }

  // Normalize LL by number of events
  double N = 1.0; //static_cast<double>(ati.numEvents(0));  // get number of events for this bin
  cout << "Number of events is: " << N << endl;
  vector<double> normLLs;
  for (const auto& p : uniqueLLs) {
    normLLs.push_back(p.second / N);
  }

  // Compute softmax using normalized LLs
  double minNormLL = *min_element(normLLs.begin(), normLLs.end());
  double maxNormLL = *max_element(normLLs.begin(), normLLs.end());
  cout << "Normalized LL range: min = " << minNormLL << ", max = " << maxNormLL << endl;

  vector<double> scores;
  double sum = 0.0;
  for (size_t i = 0; i < normLLs.size(); i++) {
    double shifted = normLLs[i] - minNormLL;
    double score = exp(-shifted / T);
    scores.push_back(score);
    sum += score;
  }
  for (auto& s : scores) s /= sum;

  // Write to file
  ofstream fout("fit_rankings.txt");
  fout << "# Temperature = " << T << "\n";
  fout << "# tag_number\tLL\tnormLL\tsoftmax_score\n";
  for (size_t i = 0; i < uniqueLLs.size(); i++) {
    fout << uniqueLLs[i].first << "\t" << uniqueLLs[i].second << "\t" << normLLs[i] << "\t" << scores[i] << "\n";
  }
  fout.close();

  // print best fit results
  // print best fit results
  if(minFitTag < 0) cout << "ALL FITS FAILED!" << endl;
  else {
    cout << "MINIMUM LIKELIHOOD FROM " << minFitTag << " of " << numRnd << " RANDOM PRODUCTION PARS = " << minLL << endl;
    gSystem->Exec(Form("cp %s_%d.fit %s.fit", fitName.data(), minFitTag, fitName.data()));
    if( seedfile.size() != 0 )
      gSystem->Exec(Form("cp %s_%d.txt %s.txt", seedfile.data(), minFitTag, seedfile.data()));
    // print summary of all fits
    summarizeFits(fitLLs);
  }

}



/*
void runBootstrapFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, int strategy, string seedfile, int numRnd, double maxFraction) {
  AmpToolsInterface ati( cfgInfo );
  string fitName = cfgInfo->fitName();

  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);
  if(strategy != 1) fitManager->setStrategy(strategy);

  vector< vector<string> > parRangeKeywords = cfgInfo->userKeywordArguments("parRange");

  // keep track of mean and std of bootstrap likelihood
  double meanLL, stdLL;

  // For each bootstrap iteration:
  for (int i = 0; i < numBootstrapItersIters; i++) {
    cout << endl << "###############################" << endl;
    cout << "BOOTSTRAP " << i << " OF " << numBootstrapItersIters << endl;
    cout << endl << "###############################" << endl;

    // 1) Create a new ROOTDataReaderBootstrap with a fresh seed
    //    (Example constructor depends on your ROOTDataReaderBootstrap interface)
    unsigned int currentSeed = seedStart + i;
    DataReader* newReader = new ROOTDataReaderBootstrap( inputFiles, currentSeed );
    
    // 2) Delete the old data reader (to avoid memory leaks),
    //    then replace it with the newly constructed one.
    //    Make sure your old pointer was dynamically allocated
    //    and is not shared by other reactions.

    DataReader* oldReader = ati.dataReader("myReaction");
    if(oldReader) delete oldReader;

    ati.m_dataReaderMap["myReaction"] = newReader; // direct map access
      // (If this map is private, you might need a setter function instead.)

    // 3) Update the pointer in the LikelihoodCalculator or reinitialize it
    //    *if* your code stores the data reader internally in the LikelihoodCalculator
    //    at construction. Typically, each LikelihoodCalculator has a pointer
    //    to the data it uses.  
    // 
    //    If there's a method like setDataReader(...), you can do:
    LikelihoodCalculator* likCalc = ati.likelihoodCalculator("myReaction");
    likCalc->setDataReader(newReader); // hypothetical method, if it exists

    //    If there's no method to change the pointer after construction,
    //    you may need to reconstruct the LikelihoodCalculator or do a partial re-init.

    // 4) Let AmpTools know it should recalculate amplitudes
    ati.invalidateAmps();

    // re-initialize parameters from configuration file (reset those not randomized)
    ati.reinitializePars();

    // Execute the fit
    if(useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if(hesse)
       fitManager->hesseEvaluation();

    bool fitFailed;
    if(strategy == 0){
      fitFailed = (fitManager->status() != 0); //eMatrixStatus not available for strategy #0
    }
    else{
      fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);
    }

    if( fitFailed )
      cout << "ERROR: fit failed use results with caution..." << endl;

    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;

    ati.finalizeFit(to_string(i));

    if( seedfile.size() != 0 && !fitFailed ){
      string seedfile_rand = seedfile + Form("_%d.txt", i);
      ati.fitResults()->writeSeed( seedfile_rand );
    }
  }

  // print bootstrap statistics
  // ... ...
}
*/



void runParScan(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, int strategy, string seedfile, string parScan) {
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
  if(strategy != 1) fitManager->setStrategy(strategy);

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
      string seedfile_scan = seedfile + Form("_scan_%d.txt", i);
      ati.fitResults()->writeSeed( seedfile_scan );
    }
  }
}

int main( int argc, char* argv[] ){

   // set default parameters

   bool useMinos = false;
   bool hesse = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   unsigned int randomSeed=static_cast<unsigned int>(time(NULL));
   int maxIter = 10000;
   int strategy = 1;
   int numBootstrapIters = 0;
   float maxFraction = 0.5;
   double temperture = 5;  // default temperature


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
      if (arg == "-st"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  strategy = atoi(argv[++i]); 
         } 
      if (arg == "-mf") {
          if ((i + 1 == argc) || (argv[i + 1][0] == '-')) {
              arg = "-h"; // Missing argument case
          } else {
              try {
                  float value = std::stof(argv[++i]); // Convert input to float
                  if (value >= 0.0 && value <= 1.0) {
                      maxFraction = value; // Assign if within range
                  } else {
                      arg = "-h"; // Out-of-range case
                  }
              } catch (const std::invalid_argument&) {
                  arg = "-h"; // Non-numeric input case
              }
          }
      }
      if (arg == "-b"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  numBootstrapIters = atoi(argv[++i]); 
         } 
      if (arg == "-H") hesse = true;
      if (arg == "-p"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  scanPar = argv[++i]; }
      if (arg == "-T") {
          if ((i + 1 == argc) || (argv[i + 1][0] == '-')) arg = "-h";
          else temperture = std::stod(argv[++i]);}
      if (arg == "-h"){
         cout << endl << " Usage for: " << argv[0] << endl << endl;
         cout << "   -b  <int> \t\t Perform <int> fits each with data sampled with replacement." << endl;
         cout << "   -st <int> \t\t Set MINUIT strategy to be <int> (default is 1, fastest 0, most accurate 2)." << endl;
         cout << "   -n \t\t\t Use MINOS instead of MIGRAD" << endl;
         cout << "   -H \t\t\t Evaluate HESSE matrix after minimization" << endl;
         cout << "   -c <file>\t\t Use config file" << endl;
         cout << "   -s <output file>\t for seeding next fit based on this fit (optional)" << endl;
         cout << "   -r <int>\t\t Perform <int> fits each seeded with random parameters" << endl;
         cout << "   -rs <int>\t\t Sets the random seed used by the random number generator for the fits with randomized initial parameters. If not set will use the time()" << endl;
         cout << "   -p <parameter> \t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
         cout << "   -m <int>\t\t Maximum number of fit iterations" << endl; 
         cout << "   -mf <float>\t\t Set maximum coherent sum fraction per amplitude (range: 0.0 to 1.0)." << endl;
         cout << "   -T <double>\t\t Softmax temperature (default = 0.1, smaller = sharper selection)" << endl;
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
   AmpToolsInterface::registerAmplitude( TwoPSMoment() );
   AmpToolsInterface::registerAmplitude( VecPSMoment() );
   AmpToolsInterface::registerAmplitude( KStarHyperon() );


   AmpToolsInterface::registerDataReader( ROOTDataReader() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderHist() );
   AmpToolsInterface::registerDataReader( FSRootDataReader() );
   AmpToolsInterface::registerDataReader( FSRootDataReaderBootstrap() );

   if (numBootstrapIters==0){
      if(numRnd==0){
          if(scanPar=="")
            runSingleFit(cfgInfo, useMinos, hesse, maxIter, strategy, seedfile);
          else
            runParScan(cfgInfo, useMinos, hesse, maxIter, strategy, seedfile, scanPar);
      } else {
          cout << "Running " << numRnd << " fits with randomized parameters with seed=" << randomSeed << endl;
          AmpToolsInterface::setRandomSeed(randomSeed);
          runRndFits(cfgInfo, useMinos, hesse, maxIter, strategy, seedfile, numRnd, maxFraction, temperture);
      }
   }else{
          cout << "Running " << numBootstrapIters << " fits with data sampled (with replacement) using seed=" << randomSeed << endl;
          cout << "This module is currently under construction." << endl;
          /*for(int i=0; i numBootstrapIters; i++) {
            cout << endl << "###############################" << endl;
            cout << "FIT " << i << " OF " << numBootstrapIters << endl;
            cout << endl << "###############################" << endl;
            ReactionInfo* rctInfo = cfgInfo->reaction(reaction);
            rctInfo->setData (classname, dataargs);
            runBootstrapFits(cfgInfo, useMinos, hesse, maxIter, strategy, seedfile, numBootstrapIters);
          }*/
   }


  return 0;
}


