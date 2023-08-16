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
    		cout << "\tomegapiDelta_ampPlotter <results file name> -o <output file name>" << endl << endl;
    		return 0;
  	}

//  	bool makePlots = true;
	string outName = "omegapiDelta_ampPlot.root";
  	string resultsName(argv[1]);
  	for (int i = 2; i < argc; i++){

    		string arg(argv[i]);

    		if (arg == "-o"){
     			outName = argv[++i];
    		}
//    		if (arg == "-n"){
//      			makePlots = false;
//    		}
    		if (arg == "-h"){
      			cout << endl << " Usage for: " << argv[0] << endl << endl;
     			cout << "\t -o <file>\t output file path" << endl;
//      			cout << "\t -n \t don't make plots "<< endl;
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

	string mJTag[] = { "p", "0", "m" };
	int mJNum = 3;

	string lTag[] = { "s", "d" };
	int lNum = 2;

	string helDelTag[] = { "m3", "m1", "p1", "p3" };
	int helDelNum = 4;

	for( int refl = 0; refl < reflNum; ++refl ){
		for( int jp = 0; jp < jpNum; ++jp ){
			for( int mJ = 0; mJ < mJNum; ++mJ ){
				for( int l = 0; l < lNum; ++l ){
					for( int helDel = 0; helDel < helDelNum; helDel++ ){
						amphistname.push_back( reflTag[refl] + jpTag[jp] + mJTag[mJ] + lTag[l] + helDelTag[helDel] );
					}
				}
			}
		}
	}

//	if( makePlots ) {

	// ************************
	// set up the plot generator
 	// ************************
	atiSetup();

	omegapi_PlotGen plotGen( results , PlotGenerator::kNoGenMC ); // slow step to load ROOT trees
	cout << " Initialized ati and PlotGen" << endl;


	TFile* plotfile = new TFile( outName.c_str(), "recreate" );
	TH1::AddDirectory( kFALSE );

	string reactionName = results.reactionList()[0];
	plotGen.enableReaction( reactionName );
	vector< string > sums = plotGen.uniqueSums();
	vector< string > amps = plotGen.uniqueAmplitudes();
	cout << "Reaction " << reactionName << " enabled with " << sums.size() << " sums and " << amps.size() << " amplitudes" << endl;

	// loop over sum configurations (one for each of the individual contributions, and the combined sum of all)

	// loop over desired amplitudes 
	for( unsigned int iamp = 0; iamp <= amphistname.size(); iamp++ ) {

		// turn all amplitudes on by default 
		for( unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {
			plotGen.enableAmp( jamp );
		}

		// turn off unwanted amplitudes and sums
      		if( iamp < amphistname.size() ) {
	
        		string locampname = amphistname[iamp];

			// turn on all sums
			for( unsigned int i = 0; i < sums.size(); i++ ) plotGen.enableSum( i );
			// turn off unwanted amplitudes
       			for( unsigned int jamp = 0; jamp < amps.size(); jamp++ ) {
				if( amps[jamp].find( locampname.data() ) == std::string::npos ) {
					plotGen.disableAmp( jamp );
	  			}
			}
		}

		unsigned int iplot = PlotGenerator::kAccMC;
		// loop over different variables
		for( unsigned int ivar = 0; ivar < OmegaPiPlotGenerator::kNumHists; ivar++ )
 
			// set unique histogram name for each plot (could put in directories...)
			string histname = plotGen.getHistogram( ivar )->name();
			if( iamp < amphistname.size() ) {
				// get name of amp for naming histogram  
				histname += "_";
				histname += amphistname[iamp];
			}

			Histogram* hist = plotGen.projection( ivar, reactionName, iplot );
			TH1* thist = hist->toRoot();
			thist->SetName( histname.c_str() );
			plotfile->cd();
			thist->Write();

		}
	}
	plotfile->Close();    
	return 0;
}

