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

#include "AMPTOOLS_DATAIO/VecPsPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/Vec_ps_moment.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

typedef VecPsPlotGenerator vecps_PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
  AmpToolsInterface::registerAmplitude( PhaseOffset() );
  AmpToolsInterface::registerAmplitude( ComplexCoeff() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( OmegaDalitz() );
  AmpToolsInterface::registerAmplitude( Vec_ps_moment() );

  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );
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
    cout << "\tvecps_plotter <results file name> -o <output file name>" << endl << endl;
    return 0;
  }

  bool showGui = false;
  string outName = "vecps_plot.root";
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

  vecps_PlotGen plotGen( results , PlotGenerator::kNoGenMC );
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

  vector<string> amphistname = {"1pps", "1p0s", "1pms", "1ppd", "1p0d", "1pmd", "1mpp", "1m0p", "1mmp", "1p", "1m"};
  vector<string> reflname = {"PosRefl", "NegRefl"};

  // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
  for (unsigned int irefl = 0; irefl <= reflname.size(); irefl++){

    // loop over desired amplitudes 
    for (unsigned int iamp = 0; iamp <= amphistname.size(); iamp++ ) {

      // turn all ampltiudes by default 
      for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {
	plotGen.enableAmp( jamp );
      }

      // turn off unwanted amplitudes and sums
      if (iamp < amphistname.size()) {
	
        string locampname = amphistname[iamp];
	
	// turn on all sums by default
	for (unsigned int i = 0; i < sums.size(); i++) plotGen.enableSum(i);

	// turn off unwanted sums for reflectivity (based on naturality)
	//cout<<"refl = "<<irefl<<endl;
	if (irefl < 2) {
	  for (unsigned int i = 0; i < sums.size(); i++){
	    if( reflname[irefl] == "NegRefl" ){
	      //ImagNegSign & RealPosSign are defined to be the negative reflectivity
	      //So, we turn off the positive reflectivity here
	      if(sums[i].find("RealNegSign") != std::string::npos || sums[i].find("ImagPosSign") != std::string::npos){
		plotGen.disableSum(i);
		//cout<<"disable sum "<<sums[i]<<"\n";
	      }
	    }
	    if( reflname[irefl] == "PosRefl" ){
	      //And, we turn off the negative reflectivity here
	      if(sums[i].find("ImagNegSign") != std::string::npos || sums[i].find("RealPosSign") != std::string::npos) {
		//cout<<"disable sum "<<sums[i]<<"\n";
		plotGen.disableSum(i);
	      }
	    }
	  }	 
	}

	// turn off unwanted amplitudes
        for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {

	  if( amps[jamp].find(locampname.data()) == std::string::npos ) {
	    plotGen.disableAmp( jamp );
	    //cout<<"disable amplitude "<<amps[jamp]<<endl;
	  }
	}

      }

      cout << "Looping over input data" << endl;
      // loop over data, accMC, and genMC
      for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
	if (iplot == PlotGenerator::kGenMC || iplot == PlotGenerator::kBkgnd) continue;
	bool singleData =  irefl == reflname.size() && iamp == amphistname.size();
	if ( iplot == PlotGenerator::kData && !singleData ) continue; // only plot data once
	
	// loop over different variables
	for (unsigned int ivar  = 0; ivar  < VecPsPlotGenerator::kNumHists; ivar++){
	  
	  // set unique histogram name for each plot (could put in directories...)
	  string histname =  "";
	  if (ivar == VecPsPlotGenerator::kVecPsMass)  histname += "MVecPs";
	  else if (ivar == VecPsPlotGenerator::kCosTheta)  histname += "CosTheta";
	  else if (ivar == VecPsPlotGenerator::kPhi)  histname += "Phi";
	  else if (ivar == VecPsPlotGenerator::kCosThetaH)  histname += "CosTheta_H";
	  else if (ivar == VecPsPlotGenerator::kPhiH)  histname += "Phi_H";
	  else if (ivar == VecPsPlotGenerator::kProd_Ang)  histname += "Prod_Ang";
	  else if (ivar == VecPsPlotGenerator::kt)  histname += "t";
	  else if (ivar == VecPsPlotGenerator::kRecoilMass)  histname += "MRecoil";
	  else if (ivar == VecPsPlotGenerator::kProtonPsMass)  histname += "MProtonPs";
	  else if (ivar == VecPsPlotGenerator::kRecoilPsMass)  histname += "MRecoilPs";
	  else if (ivar == VecPsPlotGenerator::kLambda)  histname += "Lambda";
	  else if (ivar == VecPsPlotGenerator::kDalitz)  histname += "Dalitz";
    else if (ivar == VecPsPlotGenerator::kPsiVsCosTheta) histname += "PsiVsCosTheta";
    else if (ivar == VecPsPlotGenerator::kPsiVsCosThetaH) histname += "PsiVsCosTheta_H";
    else if (ivar == VecPsPlotGenerator::kPsiVsPhiH) histname += "PsiVsPhi_H";
    else if (ivar == VecPsPlotGenerator::kProd_AngVsPhi) histname += "Prod_AngVsPhi";
    else if (ivar == VecPsPlotGenerator::kPhiVsCosTheta) histname += "PhiVsCosTheta";
    else if (ivar == VecPsPlotGenerator::kPhiHVsCosThetaH) histname += "Phi_HVsCosTheta_H";
    else if (ivar == VecPsPlotGenerator::kProtonPsMassVsCosTheta) histname += "MProtonPsVsCosTheta";
    else if (ivar == VecPsPlotGenerator::kVecPsMassVsCosTheta) histname += "MVecPsVsCosTheta";
	  else continue;	  

	  if (iplot == PlotGenerator::kData) histname += "dat";
	  if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";
	  if (iplot == PlotGenerator::kAccMC) histname += "acc";
	  if (iplot == PlotGenerator::kGenMC) histname += "gen";
	  
	  if (irefl < reflname.size()){
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
	  
	  Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
	  TH1* thist = hist->toRoot();
	  thist->SetName(histname.c_str());
	  plotfile->cd();
	  thist->Write();
	  
	}
      }
    }
  }

  plotfile->Close();

  // model parameters
  cout << "Checking Parameters" << endl;
  
  // parameters to check
  vector< string > pars;
  
  pars.push_back("dsratio");

  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "vecps_fitPars.txt" );

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

  const int nAmps = amphistname.size();
  vector<string> ampsumPosRefl[nAmps];
  vector<string> ampsumNegRefl[nAmps];
  vector< pair<string,string> > phaseDiffNames;
  
  for(unsigned int i = 0; i < fullamps.size(); i++){

    // combine amplitudes with names defined above
    for(int iamp=0; iamp<nAmps; iamp++) {
	    string locampname = amphistname[iamp];
	    
	    if(fullamps[i].find("::" + locampname) == std::string::npos) continue;
	    //cout<<locampname.data()<<" "<<fullamps[i].data()<<endl;

	    // select reflectivity
	    if(fullamps[i].find("ImagNegSign") != std::string::npos || fullamps[i].find("RealPosSign") != std::string::npos) {
	      ampsumNegRefl[iamp].push_back(fullamps[i]);
	    }
	    else {
	      ampsumPosRefl[iamp].push_back(fullamps[i]);
	    }
    }

    // second loop over amplitudes to get phase difference names
    for(unsigned int j = i+1; j < fullamps.size(); j++){

      // only keep amplitudes from the same coherent sum (and ignore constrained Real)
      if(fullamps[i].find("Real") != std::string::npos) continue;
      if(fullamps[i].find("ImagNegSign") != std::string::npos && fullamps[j].find("ImagNegSign") == std::string::npos) continue;
      if(fullamps[i].find("ImagPosSign") != std::string::npos && fullamps[j].find("ImagPosSign") == std::string::npos) continue;
	    
      phaseDiffNames.push_back( std::make_pair(fullamps[i], fullamps[j]) );

    }
  }

  for(int i = 0; i < nAmps; i++){
    if(ampsumPosRefl[i].empty()) continue;
    outfile << "FIT FRACTION (coherent sum) PosRefl " << amphistname[i] << " = "
          << results.intensity(ampsumPosRefl[i]).first / results.intensity().first << " +- "
          << results.intensity(ampsumPosRefl[i]).second / results.intensity().first << endl;
     outfile << "FIT FRACTION (coherent sum) NegRefl " << amphistname[i] << " = "
          << results.intensity(ampsumNegRefl[i]).first / results.intensity().first << " +- "
          << results.intensity(ampsumNegRefl[i]).second / results.intensity().first << endl;
  }

cout<<"Computing phase differences"<<endl;
  for(unsigned int i = 0; i < phaseDiffNames.size(); i++) {
	  pair <double, double> phaseDiff = results.phaseDiff( phaseDiffNames[i].first, phaseDiffNames[i].second );
	  outfile << "PHASE DIFF " << phaseDiffNames[i].first << " " << phaseDiffNames[i].second << " " << phaseDiff.first << " " << phaseDiff.second << endl;

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
