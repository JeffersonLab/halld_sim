#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TClass.h"
#include "TFile.h"
#include "TGClient.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotFactory.h"
#include "AmpPlotter/PlotterMainWindow.h"

#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/TwoPSMoment.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/KStarHyperon.h"

#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/TwoPsPlotGenerator.h"

#include "MinuitInterface/MinuitMinimizationManager.h"

typedef TwoPsPlotGenerator TwoPs_PlotGen;

void atiSetup(){
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( PhaseOffset() );
  AmpToolsInterface::registerAmplitude( ComplexCoeff() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( TwoPSMoment() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles() );
  AmpToolsInterface::registerAmplitude( KStarHyperon() );


  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
}

using namespace std;

void printUsage(const string& programName) {
    cout << endl
         << " *** Viewing Results Using AmpPlotter and Writing ROOT Histograms *** "
         << endl
         << endl;
    cout << "Usage:" << endl
         << endl;
    cout << "\t" << programName << " <results file name> [options]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "\t-o <file>\t Specify the output file path (default: twops_plot.root)" << endl;
    cout << "\t-g       \t Enable GUI display" << endl;
    cout << "\t-n       \t Skip plotting to save time" << endl;
    cout << "\t-h       \t Show this help message" << endl
         << endl;
}

int main(int argc, char *argv[]) {

     // ************************
    // usage
    // ************************

    if (argc < 2) {
        printUsage(argv[0] ? argv[0] : "twops_plotter");
        return 0;
    }

    // Check for "-h" or invalid argv[1]
    string arg1(argv[1]);
    if (arg1 == "-h") {
        printUsage(argv[0] ? argv[0] : "twops_plotter");
        return 0;
    }

    // Default values
    bool showGui = false;
    bool plotFlag = true;
    string outName = "twops_plot.root";
    string resultsName = argv[1]; // Assign the first argument as results file name

    // Process additional command-line arguments
    for (int i = 2; i < argc; ++i) {
        string arg(argv[i]);
        if (arg == "-g") {
            showGui = true;
        } else if (arg == "-n") {
            plotFlag = false;
        } else if (arg == "-o") {
            if (i + 1 < argc) {
                outName = argv[++i];
            } else {
                cerr << "Error: Missing argument for option '-o'" << endl;
                return 1;
            }
        } else if (arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else {
            cerr << "Error: Unknown option '" << arg << "'" << endl;
            printUsage(argv[0]);
            return 1;
        }
    }


  // ************************
  // parse the command line parameters
  // ************************

  cout << "Fit results file name    = " << resultsName << endl;
  cout << "Output file name    = " << outName << endl << endl;
  cout << "Plot Flag: " << (plotFlag ? "Enabled" : "Disabled") << endl;
  // ************************
  // load the results and display the configuration info
  // ************************

  cout << "Loading Fit results" << endl;
  FitResults results(resultsName);
  if (!results.valid()) {

    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit(1);
  }
  cout << "Fit results loaded" << endl;
  // ************************
  // set up the plot generator
  // ************************
  cout << "before atisetup();" << endl;
  atiSetup();
  cout << "Plotgen results" << endl;

  TwoPs_PlotGen plotGen(results);
  cout << " Initialized ati and PlotGen" << endl;

  // ************************
  // set up an output ROOT file to store histograms
  // ************************

 
  vector<string> amphistname = {""}; 
  // {"S0+", "S0-", 
  //   "P-1+", "P-1-", "P0+", "P0-", "P1+", "P1-",  
  //   "D0+", "D0-", "D1+", "D1-", "D2+", "D2-", "D-1+", "D-1-", "D-2+", "D-2-",
  //   "S0",   "P-1",  "P0",  "P1", "D0",  "D1",  "D2",  "D-1",  "D-2",
  //   "S",  "P",  "D"};
  vector<string> reflname = {"PosRefl", "NegRefl"};

  // string reactionName = results.reactionList()[0];
  // plotGen.enableReaction(reactionName);
  // vector<string> sums = plotGen.uniqueSums();
  // vector<string> amps = plotGen.uniqueAmplitudes();
  // cout << "Reaction " << reactionName << " enabled with " << sums.size()
  //      << " sums and " << amps.size() << " amplitudes" << endl;

  // Plot diagnostic histograms in a .root file
  if(plotFlag){
    TFile *plotfile = new TFile(outName.c_str(), "recreate");
    TH1::AddDirectory(kFALSE);

    for (const auto& reactionName : results.reactionList()) {
      plotGen.enableReaction( reactionName );
          
      cout << "Reaction " << reactionName << " enabled " <<endl;
      
      vector<string> sums = plotGen.uniqueSums();
      vector<string> amps = plotGen.uniqueAmplitudes();
      cout << "with " << sums.size() << " sums and " << amps.size() << " amplitudes" << endl; 

      // loop over sum configurations (one for each of the individual contributions,
      // and the combined sum of all)
      for (unsigned int irefl = 0; irefl <= reflname.size(); irefl++) {

        // loop over desired amplitudes
        for (unsigned int iamp = 0; iamp <= amphistname.size(); iamp++) {

          // turn all ampltiudes by default
          for (unsigned int jamp = 0; jamp < amps.size(); jamp++) {
            plotGen.enableAmp(jamp);
          }

          // turn off unwanted amplitudes and sums
          if (iamp < amphistname.size()) {

            string locampname = amphistname[iamp];

            // turn on all sums by default
            for (unsigned int i = 0; i < sums.size(); i++)
              plotGen.enableSum(i);

            // turn off unwanted sums for reflectivity (based on naturality)
            // cout<<"refl = "<<irefl<<endl;
            if (irefl < 2) {
              for (unsigned int i = 0; i < sums.size(); i++) {
                if (reflname[irefl] == "NegRefl") {
                  // ImagNegSign & RealPosSign are defined to be the negative
                  // reflectivity So, we turn off the positive reflectivity here
                  if (sums[i].find("ReZ_1+P") != std::string::npos ||
                      sums[i].find("ImZ_1-P") != std::string::npos) {
                    plotGen.disableSum(i);
                    // cout<<reflname[irefl]<<" disable sum "<<sums[i]<<"\n";
                  }
                }
                if (reflname[irefl] == "PosRefl") {
                  // And, we turn off the negative reflectivity here
                  if (sums[i].find("ImZ_1+P") != std::string::npos ||
                      sums[i].find("ReZ_1-P") != std::string::npos) {
                    plotGen.disableSum(i);
                    // cout<<reflname[irefl]<<" disable sum "<<sums[i]<<"\n";
                  }
                }
              }
            }

            // turn off unwanted amplitudes
            for (unsigned int jamp = 0; jamp < amps.size(); jamp++) {

              if (amps[jamp].find(locampname.data()) == std::string::npos) {
                plotGen.disableAmp(jamp);
                // cout<<"disable amplitude "<<amps[jamp]<<endl;
              }
            }
          }

          cout << "Looping over input data" << endl;
          // loop over data, accMC, and genMC
          for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++) {
            //if (iplot == PlotGenerator::kGenMC || iplot == PlotGenerator::kBkgnd)
            // if (iplot == PlotGenerator::kGenMC)
            //   continue;
            bool singleData =
                irefl == reflname.size() && iamp == amphistname.size();
            if (iplot == PlotGenerator::kData && !singleData)
              continue; // only plot data once

            // loop over different variables
            for (unsigned int ivar = 0; ivar < TwoPsPlotGenerator::kNumHists;
                ivar++) {

              // set unique histogram name for each plot (could put in
              // directories...)
              string histname = "";
              if (ivar == TwoPsPlotGenerator::k2PsMass)
                histname += "M2Ps";
              else if (ivar == TwoPsPlotGenerator::kLambdaKMass)
                histname += "MLambdaK";
              else if (ivar == TwoPsPlotGenerator::kLambdaPiMass)
                histname += "MLambdaPi";
              else if (ivar == TwoPsPlotGenerator::kLambdaMass)
                histname += "MLambda";
              else if (ivar == TwoPsPlotGenerator::kdaughter1Mass)
                histname += "Mdaughter1"; // Proton mass
              else if (ivar == TwoPsPlotGenerator::kdaughter2Mass)
                histname += "Mdaughter2"; // Pi- mass
              else if (ivar == TwoPsPlotGenerator::kPiCosTheta)
                histname += "CosTheta";
              else if (ivar == TwoPsPlotGenerator::kPhiK)
                histname += "PhiK";
              else if (ivar == TwoPsPlotGenerator::kPhiPi)
                histname += "PhiPi";
              else if (ivar == TwoPsPlotGenerator::kPhiLambda)
                histname += "PhiLambda";
              else if (ivar == TwoPsPlotGenerator::kThetaK)
                histname += "ThetaK";
              else if (ivar == TwoPsPlotGenerator::kThetaPi)
                histname += "ThetaPi";
              else if (ivar == TwoPsPlotGenerator::kThetaLambda)
                histname += "ThetaLambda";
              else if (ivar == TwoPsPlotGenerator::kMomK)
                histname += "MomK";
              else if (ivar == TwoPsPlotGenerator::kMomPi)
                histname += "MomPi";
              else if (ivar == TwoPsPlotGenerator::kMomLambda)
                histname += "MomLambda";
              else if (ivar == TwoPsPlotGenerator::kPhi_LAB)
                histname += "Phi_LAB";
              else if (ivar == TwoPsPlotGenerator::kPhi)
                histname += "Phi";
              else if (ivar == TwoPsPlotGenerator::kphi)
                histname += "phi";
              else if (ivar == TwoPsPlotGenerator::kPsi)
                histname += "psi";
              else if (ivar == TwoPsPlotGenerator::kt)
                histname += "t";
              else if (ivar == TwoPsPlotGenerator::kCosTheta_LambdaHel)
                histname += "cosTheta_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kphi_LambdaHel)
                histname += "phi_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kPhi_LambdaHel)
                histname += "Phi_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kPsi_LambdaHel)
                histname += "psi_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kCosThetaX_LambdaHel)
                histname += "cosThetaX_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kCosThetaY_LambdaHel)
                histname += "cosThetaY_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kCosThetaZ_LambdaHel)
                histname += "cosThetaZ_LambdaHel";
              else if (ivar == TwoPsPlotGenerator::kCosThetaX_Lambda)
                histname += "cosThetaX_Lambda";
              else if (ivar == TwoPsPlotGenerator::kCosThetaY_Lambda)
                histname += "cosThetaY_Lambda";
              else if (ivar == TwoPsPlotGenerator::kCosThetaZ_Lambda)
                histname += "cosThetaZ_Lambda";
              else
                continue;

              if (iplot == PlotGenerator::kData)
                histname += "dat";
              if (iplot == PlotGenerator::kBkgnd)
                histname += "bkgnd";
              if (iplot == PlotGenerator::kAccMC)
                histname += "acc";
              if (iplot == PlotGenerator::kGenMC)
                histname += "gen";

              if (irefl < reflname.size()) {
                // get name of sum for naming histogram
                string sumName = reflname[irefl];
                histname += "_";
                histname += sumName;
              }
              if (iamp < amphistname.size()) {
                // get name of amp for naming histogram
                histname += "_";
                histname += amphistname[iamp];
              }

              // add reaction name
              histname += "_";
              histname += reactionName;

              Histogram *hist = plotGen.projection(ivar, reactionName, iplot);
              TH1 *thist = hist->toRoot();
              thist->SetName(histname.c_str());
              plotfile->cd();
              thist->Write();
            }
          } // end of loop over plots
        } // end of loop over amplitudes
      }// end of loop over reflectivity 
    } // end of reaction loop
    plotfile->Close();
  } // end of plotFlag check


  // model parameters
  cout << "Checking Parameters" << endl;

  // parameters to check
  vector<string> pars = {
    "rho000", "rho100", "rho1m10",
    "rho111", "rho001", "rho101", "rho1m11",
    "rho102", "rho1m12",
    "Sigma", "Ox", "P", "T", "Oz"
  };

  

  // file for writing parameters (later switch to putting in ROOT file)
  size_t dotPos = resultsName.find_last_of('.');
  string outFileName = ((dotPos != string::npos) ? resultsName.substr(0, dotPos) : resultsName) + "_fitPars.txt";
  ofstream outfile(outFileName);

  for (unsigned int i = 0; i < pars.size(); i++) {
    double parValue = results.parValue(pars[i]);
    double parError = results.parError(pars[i]);
    outfile << parValue << "\t" << parError << "\t" << endl;
  }

  outfile << "TOTAL EVENTS = " << results.intensity().first << " +- "
          << results.intensity().second << endl;
  vector<string> fullamps = plotGen.fullAmplitudes();
  for (unsigned int i = 0; i < fullamps.size(); i++) {
    vector<string> useamp;
    useamp.push_back(fullamps[i]);
    outfile << "FIT FRACTION " << fullamps[i] << " = "
            << results.intensity(useamp).first / results.intensity().first
            << " +- "
            << results.intensity(useamp).second / results.intensity().first
            << endl;
  }

  const int nAmps = amphistname.size();
  vector<string> ampsumPosRefl[nAmps];
  vector<string> ampsumNegRefl[nAmps];
  vector<pair<string, string>> phaseDiffNames;

  for (unsigned int i = 0; i < fullamps.size(); i++) {

    // combine amplitudes with names defined above
    for (int iamp = 0; iamp < nAmps; iamp++) {
      string locampname = amphistname[iamp];

      if (fullamps[i].find("::" + locampname) == std::string::npos)
        continue;
      // cout<<locampname.data()<<" "<<fullamps[i].data()<<endl;

      // select reflectivity
      if (fullamps[i].find("ImZ_1-P") != std::string::npos ||
          fullamps[i].find("ReZ_1+P") != std::string::npos) {
        ampsumPosRefl[iamp].push_back(fullamps[i]);
      } else {
        ampsumNegRefl[iamp].push_back(fullamps[i]);
      }
    }

    // second loop over amplitudes to get phase difference names
    for (unsigned int j = i + 1; j < fullamps.size(); j++) {

      // only keep amplitudes from the same coherent sum (and ignore constrained
      // Real)
      if (fullamps[i].find("Re") != std::string::npos)
        continue;
      if (fullamps[i].find("ImZ_1-P") != std::string::npos &&
          fullamps[j].find("ImZ_1-P") == std::string::npos)
        continue;
      if (fullamps[i].find("ImZ_1+P") != std::string::npos &&
          fullamps[j].find("ImZ_1+P") == std::string::npos)
        continue;

      phaseDiffNames.push_back(std::make_pair(fullamps[i], fullamps[j]));
    }
  }

  // intensity +/- error
  // error = \sqrt{\Sum_{ij} (deriv[i] * deriv[j] * errorMatrix[i][j])}
  for (int i = 0; i < nAmps; i++) {
    if (!ampsumPosRefl[i].empty()) {
      outfile << "FIT FRACTION (coherent sum) PosRefl " << amphistname[i]
              << " = "
              << results.intensity(ampsumPosRefl[i]).first /
                     results.intensity().first
              << " +- "
              << results.intensity(ampsumPosRefl[i]).second /
                     results.intensity().first
              << endl;
    }
    if (!ampsumNegRefl[i].empty()) {
      outfile << "FIT FRACTION (coherent sum) NegRefl " << amphistname[i]
              << " = "
              << results.intensity(ampsumNegRefl[i]).first /
                     results.intensity().first
              << " +- "
              << results.intensity(ampsumNegRefl[i]).second /
                     results.intensity().first
              << endl;
    }
  }

  cout << "Computing phase differences" << endl;
  for (unsigned int i = 0; i < phaseDiffNames.size(); i++) {
    pair<double, double> phaseDiff =
        results.phaseDiff(phaseDiffNames[i].first, phaseDiffNames[i].second);
    outfile << "PHASE DIFF " << phaseDiffNames[i].first << " "
            << phaseDiffNames[i].second << " " << phaseDiff.first << " "
            << phaseDiff.second << endl;
  }

  // covariance matrix
  vector<vector<double>> covMatrix;
  covMatrix = results.errorMatrix();

  // ************************
  // start the GUI
  // ************************

  if (showGui) {

    cout << ">> Plot generator ready, starting GUI..." << endl;

    int dummy_argc = 0;
    char *dummy_argv[] = {};
    TApplication app("app", &dummy_argc, dummy_argv);

    gStyle->SetFillColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    gStyle->SetFillStyle(1001);
    gStyle->SetPalette(1);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameFillStyle(1001);

    cout << " Initialized App " << endl;
    PlotFactory factory(plotGen);
    cout << " Created Plot Factory " << endl;
    PlotterMainWindow mainFrame(gClient->GetRoot(), factory);
    cout << " Main frame created " << endl;

    app.Run();
    cout << " App running" << endl;
  }

  return 0;
}
