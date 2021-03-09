#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "AMPTOOLS_DATAIO/OmegaPiPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

typedef OmegaPiPlotGenerator omegapi_PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );

  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
}

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter and writing root histograms *** " << endl << endl;

  if (argc < 2){
    cout << "Usage:" << endl << endl;
    cout << "\tomegapi_plotter <results file name> -o <output file name>" << endl << endl;
    return 0;
  }

  bool showGui = false;
  string outName = "omegapi_plot.root";
  string resultsName(argv[1]);
  for (int i = 2; i < argc; i++){

    string arg(argv[i]);

    if (arg == "-g"){
      showGui = true;
    }
    if (arg == "-o"){
      outName = argv[++i];
    }
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -o <file>\t output file path" << endl;
      cout << "\t -g <file>\t show GUI" << endl;
      exit(1);
    }
  }


    // ************************
    // parse the command line parameters
    // ************************

  cout << "Fit results file name    = " << resultsName << endl;
  cout << "Output file name    = " << outName << endl << endl;

    // ************************
    // load the results and display the configuration info
    // ************************

  cout << "Loading Fit results" << endl;
  FitResults results( resultsName );
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }
   cout << "Fit results loaded" << endl;
    // ************************
    // set up the plot generator
    // ************************
	cout << "before atisetup();"<< endl;
  atiSetup();
        cout << "Plotgen results"<< endl;

  omegapi_PlotGen plotGen( results , PlotGenerator::kNoGenMC );
  cout << " Initialized ati and PlotGen" << endl;

    // ************************
    // set up an output ROOT file to store histograms
    // ************************

  TFile* plotfile = new TFile( outName.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  string reactionName = results.reactionList()[0];
  plotGen.enableReaction( reactionName );
  vector<string> sums = plotGen.uniqueSums();
  vector<string> amps = plotGen.uniqueAmplitudes();
  cout << "Reaction " << reactionName << " enabled with " << sums.size() << " sums and " << amps.size() << " amplitudes" << endl;

  const int nAmpHist = 9;
  string amphistname[nAmpHist] = {"1pps", "1p0s", "1pms", "1ppd", "1p0d", "1pmd", "1p", "1m"};

  // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
  for (unsigned int isum = 0; isum <= sums.size(); isum++){

    // turn on all sums by default
    for (unsigned int i = 0; i < sums.size(); i++){
      plotGen.enableSum(i);
    }

    // for individual contributions turn off all sums but the one of interest
    if (isum < sums.size()){
      for (unsigned int i = 0; i < sums.size(); i++){
        if (i != isum) plotGen.disableSum(i);
      }
    }

   cout << "Looping over input data" << endl;
    // loop over data, accMC, and genMC
    for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
      if (iplot == PlotGenerator::kGenMC) continue;
      if (isum < sums.size() && iplot == PlotGenerator::kData) continue; // only plot data once

      // loop over different variables
      for (unsigned int ivar  = 0; ivar  < OmegaPiPlotGenerator::kNumHists; ivar++){

      // set unique histogram name for each plot (could put in directories...)
        string histname =  "";
             if (ivar == OmegaPiPlotGenerator::kOmegaPiMass)  histname += "MOmegaPi";
             else if (ivar == OmegaPiPlotGenerator::kCosTheta)  histname += "CosTheta";
             else if (ivar == OmegaPiPlotGenerator::kPhi)  histname += "Phi";
             else if (ivar == OmegaPiPlotGenerator::kCosThetaH)  histname += "CosTheta_H";
             else if (ivar == OmegaPiPlotGenerator::kPhiH)  histname += "Phi_H";
             else if (ivar == OmegaPiPlotGenerator::kProd_Ang)  histname += "Prod_Ang";
             else if (ivar == OmegaPiPlotGenerator::kt)  histname += "t";
	     else if (ivar == OmegaPiPlotGenerator::kRecoilMass)  histname += "MRecoil";

	else continue;

        if (iplot == PlotGenerator::kData) histname += "dat";
        if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";
        if (iplot == PlotGenerator::kAccMC) histname += "acc";
        if (iplot == PlotGenerator::kGenMC) histname += "gen";

        if (isum < sums.size()){
          //ostringstream sdig;  sdig << (isum + 1);
          //histname += sdig.str();

	  // get name of sum for naming histogram
          string sumName = sums[isum];
          histname += "_";
          histname += sumName;
        }

        Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
        TH1* thist = hist->toRoot();
        thist->SetName(histname.c_str());
        plotfile->cd();
        thist->Write();

      }
    }
  }

  plotfile->Close();

  // model parameters
  cout << "Checking Parameters" << endl;
  
  // parameters to check
  vector< string > pars;
  
  pars.push_back("dalitz_alpha");
  pars.push_back("dalitz_beta");
  //pars.push_back("dalitz_gamma");
  //pars.push_back("dalitz_delta");
  pars.push_back("dsratio");

  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "omegapi_fitPars.txt" );

  for(unsigned int i = 0; i<pars.size(); i++) {
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    outfile << parValue << "\t" << parError << "\t" << endl;
  }

  outfile << "TOTAL EVENTS = " << results.intensity().first << " +- " << results.intensity().second << endl;
  vector<string> fullamps = plotGen.fullAmplitudes();
  for (unsigned int i = 0; i < fullamps.size(); i++){
    vector<string> useamp;  useamp.push_back(fullamps[i]);
    outfile << "FIT FRACTION " << fullamps[i] << " = "
         << results.intensity(useamp).first /
            results.intensity().first <<  " +- "
         << results.intensity(useamp).second /
            results.intensity().first <<  endl;
  }

  const int nAmps = 9;
  string ampname[nAmps] = {"1pps", "1p0s", "1pms", "1ppd", "1p0d", "1pmd", "1p", "1m", "2m"};
  
  vector<string> ampsumPosRefl[nAmps];
  vector<string> ampsumNegRefl[nAmps];

  for(unsigned int i = 0; i < fullamps.size(); i++){

    // combine amplitudes with names defined above
    for(int iamp=0; iamp<nAmps; iamp++) {
	    string locampname = "::" + ampname[iamp];// "::";
	    if(fullamps[i].find("NegRefl") != std::string::npos && fullamps[i].find(locampname.data()) != std::string::npos) 
		    ampsumNegRefl[iamp].push_back(fullamps[i]);
	    if(fullamps[i].find("PosRefl") != std::string::npos && fullamps[i].find(locampname.data()) != std::string::npos)
		    ampsumPosRefl[iamp].push_back(fullamps[i]);
    }
  }

  for(int i = 0; i < nAmps; i++){
    if(ampsumPosRefl[i].empty()) continue;
    outfile << "FIT FRACTION (coherent sum) PosRefl " << ampname[i] << " = "
          << results.intensity(ampsumPosRefl[i]).first / results.intensity().first << " +- "
          << results.intensity(ampsumPosRefl[i]).second / results.intensity().first << endl;
     outfile << "FIT FRACTION (coherent sum) NegRefl " << ampname[i] << " = "
          << results.intensity(ampsumNegRefl[i]).first / results.intensity().first << " +- "
          << results.intensity(ampsumNegRefl[i]).second / results.intensity().first << endl;
  }

  // covariance matrix
  vector< vector< double > > covMatrix;
  covMatrix = results.errorMatrix();

    // ************************
    // start the GUI
    // ************************

  if(showGui) {

	  cout << ">> Plot generator ready, starting GUI..." << endl;
	  
	  int dummy_argc = 0;
	  char* dummy_argv[] = {};  
	  TApplication app( "app", &dummy_argc, dummy_argv );
	  
	  gStyle->SetFillColor(10);
	  gStyle->SetCanvasColor(10);
	  gStyle->SetPadColor(10);
	  gStyle->SetFillStyle(1001);
	  gStyle->SetPalette(1);
	  gStyle->SetFrameFillColor(10);
	  gStyle->SetFrameFillStyle(1001);
   
     cout << " Initialized App " << endl;     
	  PlotFactory factory( plotGen );	
     cout << " Created Plot Factory " << endl;
	  PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
     cout << " Main frame created " << endl;
	  
	  app.Run();
     cout << " App running" << endl;
  }
    
  return 0;

}

