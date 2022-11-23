#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <cassert>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"
#include "TFile.h"

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "AMPTOOLS_DATAIO/TwoZPiPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiWt_primakoff.h"
#include "AMPTOOLS_AMPS/TwoPiWt_sigma.h"
#include "AMPTOOLS_AMPS/TwoPiW_brokenetas.h"
#include "AMPTOOLS_AMPS/TwoPitdist.h"
#include "AMPTOOLS_AMPS/TwoPiNC_tdist.h"
#include "AMPTOOLS_AMPS/TwoPiEtas_tdist.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"


using namespace std;

// for string delimiter
vector<string> stringSplit (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}



typedef TwoZPiPlotGenerator PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( TwoPiAngles() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_primakoff() );
  AmpToolsInterface::registerAmplitude( TwoPiWt_sigma() );
  AmpToolsInterface::registerAmplitude( TwoPiW_brokenetas() );
  AmpToolsInterface::registerAmplitude( TwoPitdist() );
  AmpToolsInterface::registerAmplitude( TwoPiNC_tdist() );
  AmpToolsInterface::registerAmplitude( TwoPiEtas_tdist() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
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
    cout << "\ttwopi_plotter <results file name> -o <output file name>" << endl << endl;
    return 0;
  }

  bool showGui = false;
  string outName = "twopi_plot.root";
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

  FitResults results( resultsName );
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }

    // ************************
    // set up the plot generator
    // ************************

  atiSetup();
  PlotGen plotGen( results );

    // ************************
    // set up an output ROOT file to store histograms
    // ************************

  TFile* plotfile = new TFile( outName.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  vector <string> paramlist;
  paramlist = results.ampList("Primakoff");

  vector <string> amplist;

  for (auto parsum : paramlist) {
    // cout << " sum segment=" << parsum << endl;
    vector<string> parbreak = stringSplit (parsum, "::");

    if (parbreak[1] == "Aplus"  || parbreak[1] == "IAplus"  || parbreak[1] == "EAplus" ) { 
      amplist.push_back(parbreak[2]);
      cout << " amp =" << parbreak[2] << endl;
       }
  }

  // for (auto amp : amplist) cout << " amp=" << amp << endl;
  

  string reactionName = results.reactionList()[0];
  plotGen.enableReaction( reactionName );
  vector<string> sums = plotGen.uniqueSums();

  //for (auto isum : sums) cout << " isum=" << isum << endl;

  // Enable both Aplus and Aminus sum
  plotGen.enableSum(0);
  plotGen.enableSum(1);
  // for (unsigned int isum = 0; isum < sums.size(); isum++){
  //     plotGen.disableSum(isum);
  // }



  // loop over sum, with all amplitudes turned on
  for (unsigned int iamp = 0; iamp < amplist.size(); iamp++){

    // turn on all amplist by default
    for (unsigned int i = 0; i < amplist.size(); i++){
      plotGen.enableAmp(i);
    }
  }


    // loop over sum, accMC, and genMC and kBkgnd
    for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){

      // loop over different variables
      for (unsigned int ivar  = 0; ivar  < TwoZPiPlotGenerator::kNumHists; ivar++){

        // set unique histogram name for each plot (could put in directories...)
        string histname =  "";
        if (ivar == TwoZPiPlotGenerator::k2PiMass)  histname += "M2pi";
	else if (ivar == TwoZPiPlotGenerator::kPiPCosTheta)  histname += "cosTheta";
        else if (ivar == TwoZPiPlotGenerator::kPhi)  histname += "Phi";
        else if (ivar == TwoZPiPlotGenerator::kphi)  histname += "phi";
        else if (ivar == TwoZPiPlotGenerator::kPsi)  histname += "psi";
        else if (ivar == TwoZPiPlotGenerator::kt)  histname += "t";
        else if (ivar == TwoZPiPlotGenerator::ktheta_scat)  histname += "theta_scat";
        else continue;

        if (iplot == PlotGenerator::kData) histname += "dat";
        if (iplot == PlotGenerator::kAccMC) histname += "acc";
        if (iplot == PlotGenerator::kGenMC) histname += "gen";
        if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";

        Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
        TH1* thist = hist->toRoot();
        thist->SetName(histname.c_str());
        plotfile->cd();
        thist->Write();

      }
    }



  // loop over sum, once per each amplitude
  for (unsigned int iamp = 0; iamp < amplist.size(); iamp++){

    // turn on all amplist by default
    for (unsigned int i = 0; i < amplist.size(); i++){
      plotGen.enableAmp(i);
    }

    // for individual contributions turn off all amplist but the one of interest
    for (unsigned int i = 0; i < amplist.size(); i++){
        if (i != iamp) plotGen.disableAmp(i);
      }


    // loop over data, accMC, and genMC and kBkgnd
    for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){

      // loop over different variables
      for (unsigned int ivar  = 0; ivar  < TwoZPiPlotGenerator::kNumHists; ivar++){

        // set unique histogram name for each plot (could put in directories...)
        string histname =  "";
        if (ivar == TwoZPiPlotGenerator::k2PiMass)  histname += "M2pi";
	else if (ivar == TwoZPiPlotGenerator::kPiPCosTheta)  histname += "cosTheta";
        else if (ivar == TwoZPiPlotGenerator::kPhi)  histname += "Phi";
        else if (ivar == TwoZPiPlotGenerator::kphi)  histname += "phi";
        else if (ivar == TwoZPiPlotGenerator::kPsi)  histname += "psi";
        else if (ivar == TwoZPiPlotGenerator::kt)  histname += "t";
        else if (ivar == TwoZPiPlotGenerator::ktheta_scat)  histname += "theta_scat";
        else continue;

        if (iplot == PlotGenerator::kData) histname += "dat";
        if (iplot == PlotGenerator::kAccMC) histname += "acc";
        if (iplot == PlotGenerator::kGenMC) histname += "gen";
        if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";

	  // get name of amp for naming histogram
          string ampName = amplist[iamp];
          histname += "_";
          histname += ampName;

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
    // retrieve amplitudes for output
    // ************************

  // get parameter list

  // parameters to check
  vector< string > pars;
  /* pars.push_back("Primakoff::Aplus::g1V00_re");
  pars.push_back("Primakoff::Aplus::g1V00_im");
  pars.push_back("Primakoff::Aplus::g1V11_re");
  pars.push_back("Primakoff::Aplus::g1V11_im");
  pars.push_back("Primakoff::Aplus::g1V10_re");
  pars.push_back("Primakoff::Aplus::g1V10_im");
  pars.push_back("Primakoff::Aplus::g1V1-1_re");
  pars.push_back("Primakoff::Aplus::g1V1-1_im");*/ 

  vector <string> parlist;
  parlist = results.ampList("Primakoff");
  for(unsigned int j=0; j<parlist.size(); j++) {
    cout << " j=" << j << " parlist[j]=" << parlist[j] << " " << results.realProdParName(parlist[j]) << " " << results.imagProdParName(parlist[j]) << endl;
    if ((parlist[j].find("Aplus") != string::npos) ) {
      pars.push_back(results.realProdParName(parlist[j]));
      pars.push_back(results.imagProdParName(parlist[j]));
    }
  }

  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "twopi_fitPars.txt" );
  cout << "Opened Output File twopi_fitPars.txt" << " pars.size=" << pars.size() << endl;

  for(unsigned int i = 0; i<pars.size(); i++) {
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    int ifindg = pars[i].find("g");
    outfile << pars[i].substr(ifindg) << "\t" << parValue << "\t" << parError << "\t";
  }



  // Note: For twopi_primakoff_plotter: The following computations are nonsense for amplitudes

  // covariance matrix
  vector< vector< double > > covMatrix;
  covMatrix = results.errorMatrix();

  double SigmaN = results.parValue(pars[0]) + results.parValue(pars[0]);
  double SigmaN_err = covMatrix[0][0] + covMatrix[0][0] + 2*covMatrix[0][0];

  double SigmaD = 0.5*(1 - results.parValue(pars[0])) + results.parValue(pars[0]);
  double SigmaD_err = 0.5*0.5*covMatrix[0][0] + covMatrix[0][0] - 2*0.5*covMatrix[0][0];

  double Sigma = SigmaN/SigmaD;
  double Sigma_err = fabs(Sigma) * sqrt(SigmaN_err/SigmaN/SigmaN + SigmaD_err/SigmaD/SigmaD);

  double P = 2*results.parValue(pars[0]) - results.parValue(pars[0]);
  double P_err = sqrt(2*2*covMatrix[0][0] + covMatrix[0][0] - 2*2*covMatrix[0][0]);

  Sigma = Sigma_err = P = P_err = 0;
  outfile << "Sigma" << "\t" << Sigma << "\t" << Sigma_err << "\t";
  outfile  << "P" << "\t" << P << "\t" << P_err << "\t";

  // output covariance matrix. Output only half since A+ and A- are constrained to be the same.
  for (unsigned int j=0; j< covMatrix.size()/2; j++) {
    outfile << endl;
    for (unsigned int jj=0; jj< covMatrix.size()/2; jj++) {
      outfile.width(20);
      outfile << covMatrix[j][jj];
    }
  }


  outfile << endl;

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
	  
	  PlotFactory factory( plotGen );	
	  PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
	  
	  app.Run();
  }
    
  return 0;

}

