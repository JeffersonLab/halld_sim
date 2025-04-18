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

    AmpToolsInterface::registerDataReader( ROOTDataReader() );
    AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
    AmpToolsInterface::registerDataReader( FSRootDataReader() );
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
    bool makePlots = true;
    string outName = "omegapi_plot.root";
    string resultsName(argv[1]);
    for (int i = 2; i < argc; i++){

        string arg(argv[i]);

        // if (arg == "-g"){
        //   showGui = true;
        // }
        if (arg == "-o"){
            outName = argv[++i];
        }
        if (arg == "-n"){
            makePlots = false;
        }
        if (arg == "-h"){
            cout << endl << " Usage for: " << argv[0] << endl << endl;
            cout << "\t -o <file>\t output file path" << endl;
            // cout << "\t -g \t show GUI" << endl;
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

    vector<string> amphistname = {"1ppS", "1p0S", "1pmS", "1ppD", "1p0D", "1pmD", "1mpP", "1m0P", "1mmP", "1p", "1m"};
//    vector<string> amphistname = {"0m0p", "1pps", "1p0s", "1pms", "1ppd", "1p0d", "1pmd", "1mpp", "1m0p", "1mmp", "2mp2p", "2mpp", "2m0p", "2mmp", "2mm2p", "2mp2f", "2mpf", "2m0f", "2mmf", "2mm2f", "3mp2f", "3mpf", "3m0f", "3mmf", "3mm2f", "0m", "1p", "1m", "2m", "3m"};
    vector<string> reflname = {"PosRefl", "NegRefl"};

    if(makePlots) {

        // ************************
        // set up the plot generator
        // ************************
        atiSetup();

        omegapi_PlotGen plotGen( results , PlotGenerator::kNoGenMC ); // slow step to load ROOT trees
        cout << " Initialized ati and PlotGen" << endl;
        //  plotGen.setWeightMCByIntensity( false ); // turn off weighting by model

        // ************************
        // set up an output ROOT file to store histograms
        // ************************

        TFile* plotfile = new TFile( outName.c_str(), "recreate");
        vector< TDirectory* > plotDir;

        TH1::AddDirectory(kFALSE);

        int nReactions = results.reactionList().size();

        string reactionName;
        for( int iReac = 0; iReac < nReactions; iReac++ ){
            reactionName = results.reactionList()[iReac];
            plotGen.enableReaction( reactionName );
            // create directories in plotfile
            if( iReac == 0 ){
                for( unsigned int ivar = 0; ivar < OmegaPiPlotGenerator::kNumHists; ivar++ ){
                    string dirName = plotGen.getHistogram( ivar )->name();
                    cout << "Creating directory for " << dirName.data() << " plots" << endl;
                    TDirectory *temp = plotfile->mkdir( dirName.data() );
                    plotDir.push_back( temp );
                }
            }
            vector<string> sums = plotGen.uniqueSums();
            vector<string> amps = plotGen.uniqueAmplitudes();
            cout << "Reaction " << reactionName << " enabled with " << sums.size() << " sums and " << amps.size() << " amplitudes" << endl;

            // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
            for (unsigned int irefl = 0; irefl <= reflname.size(); irefl++){

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
                        for (unsigned int i = 0; i < sums.size(); i++) plotGen.enableSum(i);

                        // turn off unwanted sums for reflectivity
                        //cout<<"refl = "<<irefl<<endl;
                        if (irefl < 2) {
                            for (unsigned int i = 0; i < sums.size(); i++){

                                bool disableSum = false;
                                if(sums[i].find("ImagNegSign") != std::string::npos || sums[i].find("RealPosSign") != std::string::npos) {
                                    if (irefl==0) disableSum = true;
                                }
                                else {
                                    if (irefl==1) disableSum = true;
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
                        if (iplot == PlotGenerator::kGenMC) continue; // || iplot == PlotGenerator::kBkgnd) continue;
                        if (irefl < reflname.size() && iamp < amphistname.size() && iplot == PlotGenerator::kData) continue; // only plot data once

                        // loop over different variables
                        for (unsigned int ivar  = 0; ivar  < OmegaPiPlotGenerator::kNumHists; ivar++){

                            // set unique histogram name for each plot (could put in directories...)
                            string histname = plotGen.getHistogram( ivar )->name();

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
                            if( nReactions > 1 ){
                                histname += "_";
                                histname += reactionName;
                            }
                            Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
                            TH1* thist = hist->toRoot();
                            thist->SetName(histname.c_str());
                            plotDir[ivar]->cd();
                            thist->Write();
                            plotfile->cd();
                        }
                    }
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

        // file for writing parameters (later switch to putting in ROOT file)
        ofstream outfile;
        outfile.open( "omegapi_fitPars.txt" );

        for(unsigned int i = 0; i<pars.size(); i++) {
            double parValue = results.parValue( pars[i] );
            double parError = results.parError( pars[i] );
            outfile << parValue << "\t" << parError << "\t" << endl;
        }
        outfile << results.covariance( "dsratio", "dphase" ) << endl;

        outfile << "TOTAL EVENTS (ACCEPTANCE CORRECTED) = " << results.intensity().first << " +- " << results.intensity().second << endl;
        outfile << "TOTAL EVENTS (NOT ACCEPTANCE CORRECTED) = " << results.intensity( false ).first << " +- " << results.intensity( false ).second << endl;
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

