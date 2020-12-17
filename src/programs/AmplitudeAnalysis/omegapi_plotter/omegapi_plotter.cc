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
#include "AMPTOOLS_AMPS/omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/vec_ps_refl.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

typedef OmegaPiPlotGenerator omegapi_PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( vec_ps_refl() );

  AmpToolsInterface::registerDataReader( ROOTDataReader() );
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

  omegapi_PlotGen plotGen( results );
  cout << " Initialized ati and PlotGen" << endl;

    // ************************
    // set up an output ROOT file to store histograms
    // ************************

  TFile* plotfile = new TFile( outName.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  string reactionName = results.reactionList()[0];
  plotGen.enableReaction( reactionName );
  vector<string> sums = plotGen.uniqueSums();
  cout << "Reaction " << reactionName << " enabled" << endl;

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

    // ************************
    // retrieve SDME parameters for plotting and asymmetry
    // ************************
/*
  cout << "Checking Parameters" << endl;
// parameters to check
  vector< string > pars;

  pars.push_back("hel_c_0_m_1");

  pars.push_back("hel_c_1_m_1");

  pars.push_back("hel_c_1_p_0");
  pars.push_back("hel_c_1_p_2");

  pars.push_back("hel_c_2_m_1");

  pars.push_back("hel_c_2_p_2");

  pars.push_back("dalitz_beta");
  pars.push_back("dalitz_gamma");
  pars.push_back("dalitz_delta");

  pars.push_back("polFrac");



  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "omegapi_fitPars.txt" );

  for(unsigned int i = 0; i<pars.size(); i++) {
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    outfile << parValue << "\t" << parError << "\t";
  }

  // covariance matrix
  vector< vector< double > > covMatrix;
  covMatrix = results.errorMatrix();

//   double SigmaN = results.parValue(pars[3]) + results.parValue(pars[6]);
//   double SigmaN_err = covMatrix[5][5] + covMatrix[8][8] + 2*covMatrix[5][8];
// 
//   double SigmaD = 0.5*(1 - results.parValue(pars[0])) + results.parValue(pars[2]);
//   double SigmaD_err = 0.5*0.5*covMatrix[2][2] + covMatrix[4][4] - 2*0.5*covMatrix[2][4];
// 
//   double Sigma = SigmaN/SigmaD;
//   double Sigma_err = fabs(Sigma) * sqrt(SigmaN_err/SigmaN/SigmaN + SigmaD_err/SigmaD/SigmaD);
//   outfile << Sigma << "\t" << Sigma_err << "\t";
// 
//   double P = 2*results.parValue(pars[6]) - results.parValue(pars[4]);
//   double P_err = sqrt(2*2*covMatrix[8][8] + covMatrix[6][6] - 2*2*covMatrix[6][8]);
//   outfile << P << "\t" << P_err << "\t";

  outfile << endl;
*/
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

