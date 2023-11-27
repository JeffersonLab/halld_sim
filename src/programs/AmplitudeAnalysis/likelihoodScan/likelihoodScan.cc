#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "TSystem.h"
#include "TH2.h"

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderWithTCut.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderHist.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_DATAIO/FSRootDataReaderTEM.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"
#include "AMPTOOLS_AMPS/DeltaAngles.h"
#include "AMPTOOLS_AMPS/SinglePS.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

void runProdParScan( ConfigurationInfo* cfgInfo, int nBins, string xAmp, double xMin, double xMax, string yAmp, double yMin, double yMax, string outfile, bool fullScan, int maxIter ) {
	AmpToolsInterface ati( cfgInfo );

	if( xMin > xMax || yMin > yMax ){
		cout << "Invalid parameter range" << endl;
		return;
	}

	TFile *fOut = TFile::Open( outfile.c_str(), "RECREATE" );

	string title = ";" + xAmp + ";" + yAmp;
	TH2F *hScan = new TH2F( "hScan", title.c_str(), nBins, xMin, xMax, nBins, yMin, yMax );

	double xBinSize = (xMax - xMin) / nBins;
	double yBinSize = (yMax - yMin) / nBins;

	string termX( "omegapi::NegHelNegPolNorm::" + xAmp );
	string termY( "omegapi::NegHelNegPolNorm::" + yAmp );

	cout << "LIKELIHOOD BEFORE MINIMIZATION = " << ati.likelihood() << endl;
	double minLL = numeric_limits< double >::max();
	
	MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
	fitManager->setMaxIterations( maxIter );

	for( int iX = 0; iX < nBins; iX++ ) {
		double realValX = xMin + ( iX + 0.5 ) * xBinSize; // use value at the center of the bin
		for( int iY = 0; iY < nBins; iY++ ) {
			if( fullScan ) ati.reinitializePars();
			ati.parameterManager()->setProductionParameter( termX, complex< double >( realValX, 0 ) );
			double realValY = yMin + ( iY + 0.5 ) * yBinSize;
			ati.parameterManager()->setProductionParameter( termY, complex< double >( realValY, 0 ) );
			bool fitFailed = false;
			if( fullScan ) {
				fitManager->migradMinimization();
				fitFailed = ( fitManager->status() != 0 || fitManager->eMatrixStatus() != 3 );
				if( fitFailed )
					cout << "ERROR: fit failed use results with caution..." << endl;
			}
			double likelihood = ati.likelihood();
			if( fullScan ) ati.finalizeFit();
			if( !fitFailed && likelihood < minLL ) minLL = likelihood;
			hScan->Fill( realValX, realValY, likelihood );
		}
	}

	cout << "MINIMUM LIKELIHOOD = " << minLL << endl;

	for( int iX = 0; iX < nBins; iX++ ) {
		double realValX = xMin + ( iX + 0.5 ) * xBinSize; // use value at the center of the bin
		for( int iY = 0; iY < nBins; iY++ ) {
			double realValY = yMin + ( iY + 0.5 ) * yBinSize;
			int bin = hScan->FindBin( realValX, realValY );
			double oldContent = hScan->GetBinContent( bin );
			double newContent = oldContent - minLL;
			hScan->SetBinContent( bin, newContent );
		}
	}


	fOut->WriteObject( hScan, "hScan" );
	fOut->Close();
}

void runParScan( ConfigurationInfo* cfgInfo, int nBins, string xAmp, double xMin, double xMax, string yAmp, double yMin, double yMax, string outfile, bool fullScan, int maxIter ) {
	AmpToolsInterface ati( cfgInfo );
	
	if( xMin > xMax || yMin > yMax ){
		cout << "Invalid parameter range" << endl;
		return;
	}

	TFile *fOut = TFile::Open( outfile.c_str(), "RECREATE" );

	string title = ";" + xAmp + ";" + yAmp;
	TH2F* hScan = new TH2F( "hScan", title.c_str(), nBins, xMin, xMax, nBins, yMin, yMax );

	double xBinSize = (xMax - xMin) / nBins;
	double yBinSize = (yMax - yMin) / nBins;

//	string termX( xAmp );
//	string termY( "omegapi::NegHelNegPolNorm::" + yAmp );

	cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;

	ParameterManager* parManager = ati.parameterManager();
	MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
	fitManager->setMaxIterations( maxIter );

	double minLL = numeric_limits<double>::max();

	for( int iX = 0; iX < nBins; iX++ ){
		double valX = xMin + ( iX + 0.5 ) * xBinSize;
		for( int iY = 0; iY < nBins; iY++ ){
			if( fullScan ) ati.reinitializePars();
			parManager->setAmpParameter( xAmp, valX );
			double valY = yMin + ( iY + 0.5 ) * yBinSize;
			parManager->setAmpParameter( yAmp, valY );

			bool fitFailed = false;
			if( fullScan ){
				fitManager->migradMinimization();
				fitFailed = ( fitManager->status() != 0 || fitManager->eMatrixStatus() != 3 );
				if( fitFailed )
					cout << "ERROR: fit failed use results with caution..." << endl;
			}
			double likelihood = ati.likelihood();
			if( fullScan ) ati.finalizeFit();
			if( !fitFailed && likelihood < minLL ) minLL = likelihood;

			hScan->Fill( valX, valY, likelihood );			
		}
	}	

	cout << "MINIMUM LIKELIHOOD = " << minLL << endl;

	for( int iX = 0; iX < nBins; iX++ ){
		double valX = xMin + ( iX + 0.5 ) * xBinSize;
		for( int iY = 0; iY < nBins; iY++ ){
			double valY = yMin + ( iY + 0.5 ) * yBinSize;
			int bin = hScan->FindBin( valX, valY );
			double oldContent = hScan->GetBinContent( bin );
			double newContent = oldContent - minLL;
			hScan->SetBinContent( bin, newContent );
		}
	}

	fOut->WriteObject( hScan, "hScan" );
	fOut->Close();
}

int main( int argc, char* argv[] ){

	// set default parameters

	string configfile;
	string outfile;
	int nBins = 10;
	bool fullScan = false;
	bool scanAmpPars = false;
	double xMin = 0.;
	double xMax = 90.;
	double yMin = 0.;
	double yMax = 90.;
	string xAmp = "m1p0sp1";
	string yAmp = "m1p0sm1";
	int maxIter = 10000;

	// parse command line

	for( int i = 1; i < argc; i++ ){

		string arg( argv[i] );

		if( arg == "-c" ){  
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else  configfile = argv[++i]; }
		if( arg == "-o" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else outfile = argv[++i]; }
		if( arg == "-b" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else nBins = atoi( argv[++i] ); }
		if( arg == "-f" ){
			fullScan = true; }
		if( arg == "-m" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else maxIter = atoi( argv[++i] ); }
		if( arg == "-x" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else xAmp = argv[++i]; }	
		if( arg == "-y" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
			else yAmp = argv[++i]; }
		if( arg == "-xMin" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
      			else xMin = atof( argv[++i] ); }
		if( arg == "-xMax" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
      			else xMax = atof( argv[++i] ); }
		if( arg == "-yMin" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
      			else yMin = atof( argv[++i] ); }
		if( arg == "-yMax" ){
			if( ( i+1 == argc ) || ( argv[i+1][0] == '-' ) ) arg = "-h";
      			else yMax = atof( argv[++i] ); }
		if( arg == "-xMinNeg" ){
			xMin *= -1.; }
		if( arg == "-xMaxNeg" ){
			xMax *= -1.; }
		if( arg == "-yMinNeg" ){
			yMin *= -1.; }
		if( arg == "-yMaxNeg" ){
			yMax *= -1.; }
		if( arg == "-ap" ){
			scanAmpPars = true; }
		if( arg == "-h" ){
         		cout << endl << " Usage for: " << argv[0] << endl << endl;
         		cout << "   -c <file>\t\t\t\t config file" << endl;
			cout << "   -o <file>\t\t\t\t output file" << endl;
			cout << "   -b <int>\t\t\t\t number of bins on each axis (default 10)" << endl;
			cout << "   -f \t\t\t\t\t run a fit at each scan point (default will simply calculate likelihood)" << endl;
			cout << "   -m <int>\t\t\t\t maximum number of fit iterations (default 10000)" << endl;
			cout << "   -x <string>\t\t\t\t amplitude on x-axis (default m1p0sp1)" << endl;
			cout << "   -y <string>\t\t\t\t amplitude on y-axis (default m1p0sm1)" << endl;
			cout << "   -xMin <double>\t\t\t minimum value on x-axis (default 0)" << endl;
			cout << "   -xMax <double>\t\t\t maximum value on x-axis (default 90)" << endl;
			cout << "   -yMin <double>\t\t\t minimum value on y-axis (default 0)" << endl;
			cout << "   -yMax <double>\t\t\t maximum value on y-axis (default 90)" << endl;
			cout << "   -xMinNeg \t\t\t\t multiply xMin by -1" << endl;
			cout << "   -xMaxNeg \t\t\t\t multiply xMax by -1" << endl;
			cout << "   -yMinNeg \t\t\t\t multiply yMin by -1" << endl;
			cout << "   -yMaxNeg \t\t\t\t multiply yMax by -1" << endl;
			cout << "   -ap \t\t\t\t\t scan amplitude parameters (default will scan production parameters)" << endl; 	
			exit(1); }
   	}

   	if (configfile.size() == 0){
      		cout << "No config file specified" << endl;
            	exit(1);
   	}

   	ConfigFileParser parser(configfile);
   	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
   	cfgInfo->display();

   	AmpToolsInterface::registerAmplitude( BreitWigner() );
   	AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
   	AmpToolsInterface::registerAmplitude( PhaseOffset() );
   	AmpToolsInterface::registerAmplitude( ComplexCoeff() );
   	AmpToolsInterface::registerAmplitude( OmegaDalitz() );
   	AmpToolsInterface::registerAmplitude( Piecewise() );
   	AmpToolsInterface::registerAmplitude( LowerVertexDelta() );
   	AmpToolsInterface::registerAmplitude( DeltaAngles() );
   	AmpToolsInterface::registerAmplitude( SinglePS() );

   	AmpToolsInterface::registerDataReader( ROOTDataReader() );
   	AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   	AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
   	AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
   	AmpToolsInterface::registerDataReader( ROOTDataReaderHist() );
   	AmpToolsInterface::registerDataReader( FSRootDataReader() );
   	AmpToolsInterface::registerDataReader( FSRootDataReaderTEM() );

	if( scanAmpPars )
		runParScan( cfgInfo, nBins, xAmp, xMin, xMax, yAmp, yMin, yMax, outfile, fullScan, maxIter );
	else	
		runProdParScan( cfgInfo, nBins, xAmp, xMin, xMax, yAmp, yMin, yMax, outfile, fullScan, maxIter );

  	return 0;
}


