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
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"

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
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( PhaseOffset() );

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
  bool makePlots = false;
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
    if (arg == "-p"){
      makePlots = true;
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

   vector<string> amphistname = {"0m0p", "1pps", "1p0s", "1pms", "1ppd", "1p0d", "1pmd", "1mpp", "1m0p", "1mmp", "2mp2p", "2mpp", "2m0p", "2mmp", "2mm2p", "2mp2f", "2mpf", "2m0f", "2mmf", "2mm2f", "3mp2f", "3mpf", "3m0f", "3mmf", "3mm2f", "0m", "1p", "1m", "2m", "3m"};
  vector<string> reflname = {"PosRefl", "NegRefl"};

if( makePlots ) {
    // ************************
    // set up the plot generator
    // ************************
	cout << "before atisetup();"<< endl;
  atiSetup();
        cout << "Plotgen results"<< endl;

	cout << "before atisetup();"<< endl;
  atiSetup();
        cout << "Plotgen results"<< endl;

  omegapi_PlotGen plotGen( results , PlotGenerator::kNoGenMC ); // slow step to load ROOT trees...
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
	
	// parse amplitude name for naturality to compute reflectivity 
	int j = locampname[0]-'0';
	int parity = 0;
	if( locampname[1] == 'p' ) parity = +1;
	else if( locampname[1] == 'm' ) parity = -1;
	else cout<<"Undefined parity in amplitude"<<endl;
	int naturality = parity*pow(-1,j);

	// turn on all sums by default
	for (unsigned int i = 0; i < sums.size(); i++) plotGen.enableSum(i);

	// turn off unwanted sums for reflectivity (based on naturality)
	//cout<<"refl = "<<irefl<<endl;
	if (irefl < 2) {
	  for (unsigned int i = 0; i < sums.size(); i++){
	    
	    bool disableSum = false;
	    if(sums[i].find("ImagNegSign") != std::string::npos || sums[i].find("RealPosSign") != std::string::npos) {
	      if (naturality>0 && irefl==1) disableSum = true;
	      if (naturality<0 && irefl==0) disableSum = true;
	    }
	    else {
	      if (naturality>0 && irefl==0) disableSum = true;
	      if (naturality<0 && irefl==1) disableSum = true;
	    }
	    
	    if(disableSum) {
	      plotGen.disableSum(i);
	      //cout<<"disable sum "<<sums[i]<<endl;
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
	if (irefl < reflname.size() && iamp < amphistname.size() && iplot == PlotGenerator::kData) continue; // only plot data once
	
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
	  else if (ivar == OmegaPiPlotGenerator::kProtonPiMass)  histname += "MProtonPi";
	  else if (ivar == OmegaPiPlotGenerator::kRecoilPiMass)  histname += "MRecoilPi";
	  
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
}

  // model parameters
  cout << "Checking Parameters" << endl;
  
  // parameters to check
  vector< string > pars;
  
  pars.push_back("dalitz_alpha");
  pars.push_back("dalitz_beta");
  //pars.push_back("dalitz_gamma");
  //pars.push_back("dalitz_delta");
  pars.push_back("dsratio");
  pars.push_back("dphase");

  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "omegapi_fitPars.txt" );

  for(unsigned int i = 0; i<pars.size(); i++) {
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    outfile << parValue << "\t" << parError << "\t" << endl;
  }

  outfile << "TOTAL EVENTS = " << results.intensity().first << " +- " << results.intensity().second << endl;
  vector<string> fullamps = results.ampList(); //plotGen.fullAmplitudes();
  for (unsigned int i = 0; i < fullamps.size(); i++){
    //cout<<fullamps[i].data()<<endl;
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

	    // parse amplitude name for naturality to compute reflectivity 
	    int j = locampname[0]-'0';
	    int parity = 0;
	    if( locampname[1] == 'p' ) parity = +1;
	    else if( locampname[1] == 'm' ) parity = -1;
	    else cout<<"Undefined parity in amplitude"<<endl;
	    int naturality = parity*pow(-1,j);	   
 
	    // select reflectivity (based on naturality)
	    if(fullamps[i].find("ImagNegSign") != std::string::npos || fullamps[i].find("RealPosSign") != std::string::npos) {
		    if (naturality>0) {
			    ampsumPosRefl[iamp].push_back(fullamps[i]);
		    }
		    if (naturality<0) {
			    ampsumNegRefl[iamp].push_back(fullamps[i]);
		    }
	    }
	    else {
		    if (naturality>0) {
			    ampsumNegRefl[iamp].push_back(fullamps[i]);
		    }
		    if (naturality<0) {
			    ampsumPosRefl[iamp].push_back(fullamps[i]);
		    }
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

  cout<<"Computing fit fractions"<<endl;
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

  ///////////////////////////////////////////////////////////////////////////////
  // collect ensemble of fit results (random starting parameters or bootstrap) //
  ///////////////////////////////////////////////////////////////////////////////

  cout << "Writing file with ensemble of ";

  // loop over ensemble of results
  int maxFits = 25;
  const int nFullAmps = fullamps.size();
  const int nPhaseDiff = phaseDiffNames.size();
  vector<double> likelihood;
  vector<double> fitFraction[nFullAmps], fitFractionError[nFullAmps];
  vector<double> fitFractionCoherentPosRefl[nAmps], fitFractionCoherentPosReflError[nFullAmps];
  vector<double> fitFractionCoherentNegRefl[nAmps], fitFractionCoherentNegReflError[nFullAmps];
  vector<double> phaseDiff[nPhaseDiff], phaseDiffError[nPhaseDiff];
  for(int i=0; i<maxFits; i++){
	  TString resultName = resultsName; resultName.ReplaceAll(".fit", Form("_%d.fit", i));
	  FitResults result(resultName.Data());
	  
	  // check if fits 
	  if( !result.valid() ) continue;
	  if( result.lastMinuitCommandStatus() != 0) continue; 

	  likelihood.push_back(result.likelihood());
	  for (unsigned int i = 0; i < fullamps.size(); i++){
		  vector<string> useamp;  useamp.push_back(fullamps[i]);
		  
		  fitFraction[i].push_back( result.intensity(useamp).first / result.intensity().first );
		  fitFractionError[i].push_back( result.intensity(useamp).second / result.intensity().first );	  
	  }

	  for(int i = 0; i < nAmps; i++){
		  if(ampsumPosRefl[i].empty()) continue;
		  fitFractionCoherentPosRefl[i].push_back(result.intensity(ampsumPosRefl[i]).first / result.intensity().first);
		  fitFractionCoherentPosReflError[i].push_back(result.intensity(ampsumPosRefl[i]).second / result.intensity().first);
		  fitFractionCoherentNegRefl[i].push_back(result.intensity(ampsumNegRefl[i]).first / result.intensity().first);
		  fitFractionCoherentNegReflError[i].push_back(result.intensity(ampsumNegRefl[i]).second / result.intensity().first);	  
	  }

	  for(unsigned int i = 0; i < phaseDiffNames.size(); i++) {
		  pair <double, double> phase = result.phaseDiff( phaseDiffNames[i].first, phaseDiffNames[i].second );
		  phaseDiff[i].push_back(phase.first);
		  phaseDiffError[i].push_back(phase.second);	  
	  }
  }

  cout << likelihood.size() << endl;

  // file for writing parameters
  ofstream outfile_ensemble;
  outfile_ensemble.open( "omegapi_ensemble.txt" );

  outfile_ensemble << "TOTAL EVENTS = " << results.intensity().first << ", " << results.intensity().second << endl;
  outfile_ensemble << "LIKELIHOOD     ";
  for(size_t j=0; j<likelihood.size(); j++) outfile_ensemble << ", " << std::setprecision(12) << likelihood[j];
  outfile_ensemble << endl;

  for (unsigned int i = 0; i < fullamps.size(); i++){
	  vector<string> useamp;  useamp.push_back(fullamps[i]);
	  outfile_ensemble << "FIT FRACTION     " << fullamps[i] << " ";
	  for(size_t j=0; j<fitFraction[i].size(); j++) outfile_ensemble << ", " << fitFraction[i][j];
	  outfile_ensemble << endl;
	  outfile_ensemble << "FIT FRACTION ERR " << fullamps[i] << " ";
	  for(size_t j=0; j<fitFractionError[i].size(); j++) outfile_ensemble << ", " << fitFractionError[i][j];
	  outfile_ensemble << endl;
  }
  
  for(int i = 0; i < nAmps; i++){
	  if(ampsumPosRefl[i].empty()) continue;
	  outfile_ensemble << "FIT FRACTION (coherent sum) PosRefl     " << amphistname[i] << " ";
	  for(size_t j=0; j<fitFractionCoherentPosRefl[i].size(); j++) outfile_ensemble << ", " << fitFractionCoherentPosRefl[i][j];
	  outfile_ensemble << endl;
	  outfile_ensemble << "FIT FRACTION (coherent sum) PosRefl ERR " << amphistname[i] << " ";
	  for(size_t j=0; j<fitFractionCoherentPosReflError[i].size(); j++) outfile_ensemble << ", " << fitFractionCoherentPosReflError[i][j];
	  outfile_ensemble << endl;
	  outfile_ensemble << "FIT FRACTION (coherent sum) NegRefl     " << amphistname[i] << " ";
	  for(size_t j=0; j<fitFractionCoherentNegRefl[i].size(); j++) outfile_ensemble << ", " << fitFractionCoherentNegRefl[i][j];
	  outfile_ensemble << endl;
	  outfile_ensemble << "FIT FRACTION (coherent sum) NegRefl ERR " << amphistname[i] << " ";
	  for(size_t j=0; j<fitFractionCoherentNegReflError[i].size(); j++) outfile_ensemble << ", " << fitFractionCoherentNegReflError[i][j];
	  outfile_ensemble <<endl;
  }

  for(unsigned int i = 0; i < phaseDiffNames.size(); i++) {
	  outfile_ensemble << "PHASE DIFF     " << phaseDiffNames[i].first << " " << phaseDiffNames[i].second << " "; 
	  for(size_t j=0; j<phaseDiff[i].size(); j++) outfile_ensemble << ", " << phaseDiff[i][j];
	  outfile_ensemble << endl;
	  outfile_ensemble << "PHASE DIFF ERR " << phaseDiffNames[i].first << " " << phaseDiffNames[i].second << " "; 
	  for(size_t j=0; j<phaseDiffError[i].size(); j++) outfile_ensemble << ", " << phaseDiffError[i][j];
	  outfile_ensemble << endl;
  }

  // covariance matrix
  vector< vector< double > > covMatrix;
  covMatrix = results.errorMatrix();

    // ************************
    // start the GUI
    // ************************
  /*
  if(makePlots && showGui) {

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
     cout << " App running" << endl;
  }
  */
    
  return 0;

}

