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
#include "AMPTOOLS_DATAIO/FSRootDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_AMPS/omegapi_amplitude.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"
#include "AMPTOOLS_AMPS/DeltaAngles.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

typedef OmegaPiPlotGenerator omegapi_PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
  AmpToolsInterface::registerAmplitude( PhaseOffset() );
  AmpToolsInterface::registerAmplitude( ComplexCoeff() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( OmegaDalitz() );
  AmpToolsInterface::registerAmplitude( LowerVertexDelta() );
  AmpToolsInterface::registerAmplitude( DeltaAngles() );

  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );
  AmpToolsInterface::registerDataReader( FSRootDataReaderTEM() );
}

using namespace std;

int main( int argc, char* argv[] ){


	// ************************
    	// usage
    	// ************************

  	cout << endl << " *** Viewing Results Using AmpPlotter and writing root histograms *** " << endl << endl;

  	if (argc < 2){
    		cout << "Usage:" << endl << endl;
    		cout << "\tomegapiDelta_plotter <results file name> -o <output file name>" << endl << endl;
    		return 0;
  	}

  	bool showGui = false;
  	bool makePlots = true;
  	string outName = "omegapiDelta_plot.root";
  	string resultsName(argv[1]);
  	for (int i = 2; i < argc; i++){

    		string arg(argv[i]);

    		if (arg == "-g"){
      			showGui = true;
    		}
    		if (arg == "-o"){
      			outName = argv[++i];
    		}
    		if (arg == "-n"){
      			makePlots = false;
    		}
    		if (arg == "-h"){
      			cout << endl << " Usage for: " << argv[0] << endl << endl;
      			cout << "\t -o <file>\t output file path" << endl;
      			cout << "\t -g \t show GUI" << endl;
      			cout << "\t -n \t don't make plots "<< endl;
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

	vector< string > amphistname;

	string reflTag[] = { "p", "m" };
	int reflNum = 2;

	string jpTag[] = { "1p" };
	int jpNum = 1;

	string mJTag[] = { "m", "0", "p" };
	int mJNum = 3;

	string lTag[] = { "s" };
	int lNum = 1;

	string helDelTag[] = { "m3", "m1", "p1", "p3" };
	int helDelNum = 4;

	for( int refl = 0; refl < reflNum; ++refl ){
		for( int jp = 0; jp < jpNum; ++jp ){
			amphistname.push_back( reflTag[refl] + jpTag[jp] );
			for( int mJ = 0; mJ < mJNum; ++mJ ){
				for( int l = 0; l < lNum; ++l ){
					amphistname.push_back( reflTag[refl] + jpTag[jp] + mJTag[mJ] + lTag[l] );
					for( int helDel = 0; helDel < helDelNum; helDel++ ){
						amphistname.push_back( reflTag[refl] + jpTag[jp] + mJTag[mJ] + lTag[l] + helDelTag[helDel] );
						string tempmJ = reflTag[refl] + jpTag[jp] + lTag[l] +helDelTag[helDel];
						if( std::find( amphistname.begin(), amphistname.end(), tempmJ ) == amphistname.end() ) {
							amphistname.push_back( tempmJ );
						}
					}
				}
			}
		}
	}

	if( makePlots ) {

		// ************************
		// set up the plot generator
 		// ************************
  		atiSetup();

  		omegapi_PlotGen plotGen( results , PlotGenerator::kNoGenMC ); // slow step to load ROOT trees
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

//  	// loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
//  	for (unsigned int irefl = 0; irefl < reflNum; irefl++){

    		// loop over desired amplitudes 
    		for (unsigned int iamp = 0; iamp <= amphistname.size(); iamp++ ) {

      			// turn all ampltiudes on by default 
      			for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {
				plotGen.enableAmp( jamp );
      			}

			// turn off unwanted amplitudes and sums
	      		if (iamp < amphistname.size()) {
	
	        		string locampname = amphistname[iamp];

				// turn on all sums by default
				for (unsigned int i = 0; i < sums.size(); i++) plotGen.enableSum( i );
/*
				// pick PosHelPosPolConj sum to stay on, others get turned off
				for (unsigned int i = 0; i < sums.size(); i++){
					bool disableSum = false;
					if( sums[i].find( "PosHelPosPolConj" ) == std::string::npos ) disableSum = true;
					if( disableSum ){
						plotGen.disableSum( i );
						cout << "Disable sum: " << sums[i] << endl;
					}
				}
*/
				// turn off unwanted amplitudes
				for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {
					if( amps[jamp].find( locampname.data() ) == std::string::npos ) {
						plotGen.disableAmp( jamp );
//						cout << "Disable amplitude " << amps[jamp] << endl;
					}
				}
			}

	      		cout << "Looping over input data" << endl;
      			// loop over data, accMC, and genMC
      			for( unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++ ){
				if( iplot == PlotGenerator::kGenMC || iplot == PlotGenerator::kBkgnd ) continue; // can turn bkgnd on again when relevant...
				if( iamp < amphistname.size() && iplot == PlotGenerator::kData ) continue; // only plot data once
	
				// loop over different variables
				for( unsigned int ivar = 0; ivar < OmegaPiPlotGenerator::kNumHists; ivar++ ){
	  
	  				// set unique histogram name for each plot (could put in directories...)
	  				string histname = plotGen.getHistogram( ivar )->name();

	  				if (iplot == PlotGenerator::kData) histname += "dat";
	  				if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";
	  				if (iplot == PlotGenerator::kAccMC) histname += "acc";
	  				if (iplot == PlotGenerator::kGenMC) histname += "gen";

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
	  	pars.push_back("dphase");

		// full amplitude names	
		vector< string > fullAmps = plotGen.fullAmplitudes();

  	// ***************************
  	// .csv file of fit parameters
  	// ***************************

	/* Note that code convention below is that every written entry starts with a comma, rather than ends with it.
	   This avoids having to "detect" if an entry is the last of the line
	*/
/*
	// create csvFile
	ofstream csvFile;
	csvFile.open( "omegapiDelta_fitPars.csv" );

	// write non-amplitude headers first
	csvFile << "eMatrixStatus" << ",lastMinuitCommandStatus" << ",likelihood" << ",total_events" << ",total_events_err";

	for( string par : pars ){
		csvFile << "," << par << "," << par + "_err";
	}

	// set up map of amplitude based header name, to a vector that stores the full amplitude name
	// Excluded variables mean they are summed over
	// ZB: sum of all JP = 1+, l = 0 waves is 1ps
  	std::map< string, vector<string> > ampSum_eJPmJlD;
  	std::map< string, vector<string> > ampSum_JPmJlD;
  	std::map< string, vector<string> > ampSum_JPmJD;
  	std::map< string, vector<string> > ampSum_JPlD;
  	std::map< string, vector<string> > ampSum_eJPmJD;
  	std::map< string, vector<string> > ampSum_eJPlD;
  	std::map< string, vector<string> > ampSum_eJPD;
  	std::map< string, vector<string> > ampSum_JPD;
  	std::map< string, vector<string> > ampSum_eJPmJl;
  	std::map< string, vector<string> > ampSum_JPmJl;
  	std::map< string, vector<string> > ampSum_JPmJ;
  	std::map< string, vector<string> > ampSum_JPl;
  	std::map< string, vector<string> > ampSum_eJPmJ;
  	std::map< string, vector<string> > ampSum_eJPl;
  	std::map< string, vector<string> > ampSum_eJP;
  	std::map< string, vector<string> > ampSum_JP;
  	std::map< string, pair<string,string> > ampPhaseDiffs; // eJPmlD_eJPmlD

	// WRITE HEADERS
	for( unsigned int i = 0; i < fullAmps.size(); i++ ){
		string fullAmp = fullAmps[i];
		string delim = "::";

		// extract eJPmJlD amplitude
		// assumes fullAmps always in "reaction::sum::amp" format
		size_t position1 = fullAmp.find( delim ) + delim.length();
		size_t positionAmp = fullAmp.find( delim, position1 ) + delim.length();
		string eJPmJlD = fullAmp.substr( positionAmp, fullAmp.length() - positionAmp );

		string e = eJPmJlD.substr( 0, 1 ), JP = eJPmJlD.substr( 1, 2 ), mJ = eJPmJlD.substr( 3, 2 ), l = eJPmJlD.substr( 5, 1 ), D = eJPmJlD.substr( 6, 2 );

		ampSum_eJPmJlD[eJPmJlD].push_back( fullAmp );
	}
*/
 	 	// file for writing parameters (later switch to putting in ROOT file)
 	 	ofstream outfile;
 	 	outfile.open( "omegapiDelta_fitPars.txt" );

	  	for(unsigned int i = 0; i<pars.size(); i++) {
	    		double parValue = results.parValue( pars[i] );
	    		double parError = results.parError( pars[i] );
	    		outfile << parValue << "\t" << parError << "\t" << endl;
	  	}
	  	outfile << results.covariance( "dsratio", "dphase" ) << endl;

	  	outfile << "TOTAL EVENTS = " << results.intensity().first << " +- " << results.intensity().second << endl;
	  	vector<string> fullamps = results.ampList();
	  	for (unsigned int i = 0; i < fullamps.size(); i++){
	    		vector<string> useamp;  useamp.push_back(fullamps[i]);
	    		outfile << "FIT FRACTION " << fullamps[i] << " = "
	         		<< results.intensity(useamp).first /
	            		   results.intensity().first <<  " +- "
	         		<< results.intensity(useamp).second /
	         		   results.intensity().first <<  endl;
	  	}

	  	const int nAmps = amphistname.size();
//  	vector<string> ampsumPosRefl[nAmps];
// 	vector<string> ampsumNegRefl[nAmps];
		vector< string > ampSum[nAmps];
	  	vector< pair< string, string > > phaseDiffNames;

 	 	for( unsigned int i = 0; i < fullamps.size(); i++ ){

 	   		// combine amplitudes with names defined above
 	   		for( int iamp = 0; iamp < nAmps; iamp++ ) {
				string locampname = amphistname[iamp];
	    
		    		if( fullamps[i].find("::" + locampname) == std::string::npos ) continue;
//	    			cout << locampname.data() << " " << fullamps[i].data() << endl;
	 			ampSum[iamp].push_back( fullamps[i] );
	    		// select reflectivity 
//	    		if(fullamps[i].find("ImagNegSign") != std::string::npos || fullamps[i].find("RealPosSign") != std::string::npos) {
//	      			ampsumNegRefl[iamp].push_back(fullamps[i]);
//          		}
//	    		else {
//	      			ampsumPosRefl[iamp].push_back(fullamps[i]);
//	    		}
    			}
    
 	   		// second loop over amplitudes to get phase difference names
 	   		for( unsigned int j = i+1; j < fullamps.size(); j++ ){

		    		// only keep amplitudes from one coherent sum (PosHelNegPolConj)
		    		if( fullamps[i].find( "PosHelNegPolConj" ) == std::string::npos ) continue;
		    		if( fullamps[j].find( "PosHelNegPolConj" ) == std::string::npos ) continue;
//	    		if(fullamps[i].find("ImagNegSign") != std::string::npos && fullamps[j].find("ImagNegSign") == std::string::npos) continue;
//	    		if(fullamps[i].find("ImagPosSign") != std::string::npos && fullamps[j].find("ImagPosSign") == std::string::npos) continue;
	    
			    	phaseDiffNames.push_back( std::make_pair( fullamps[i], fullamps[j] ) );
		    	}
  		}

  		cout << "Computing fit fractions" << endl;
  		for( int i = 0; i < nAmps; i++ ){
    			if( ampSum[i].empty()) continue;
    			outfile << "FIT FRACTION (coherent sum) " << amphistname[i] << " = "
          			<< results.intensity( ampSum[i] ).first / results.intensity().first << " +- "
          			<< results.intensity( ampSum[i] ).second / results.intensity().first << endl;
//     		outfile << "FIT FRACTION (coherent sum) NegRefl " << amphistname[i] << " = "
//          		<< results.intensity(ampsumNegRefl[i]).first / results.intensity().first << " +- "
//          		<< results.intensity(ampsumNegRefl[i]).second / results.intensity().first << endl;
  		}

 	 	cout << "Computing phase differences" << endl;
  		for( unsigned int i = 0; i < phaseDiffNames.size(); i++ ) {
		  	pair< double, double > phaseDiff = results.phaseDiff( phaseDiffNames[i].first, phaseDiffNames[i].second );
		  	outfile << "PHASE DIFF " << phaseDiffNames[i].first << " " << phaseDiffNames[i].second << " " << phaseDiff.first << " " << phaseDiff.second << endl;

  		}

  		// covariance matrix
  		vector< vector< double > > covMatrix;
  		covMatrix = results.errorMatrix();

    	// ************************
    	// start the GUI
    	// ************************

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
	  	gStyle->SetFrameFillStyle(1001);
   
     		cout << " Initialized App " << endl;     
	  	PlotFactory factory( plotGen );	
     		cout << " Created Plot Factory " << endl;
	  	PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
     		cout << " Main frame created " << endl;
	  
	  	app.Run();
	     	cout << " App running" << endl;
  	}

	}    
  	return 0;

}

