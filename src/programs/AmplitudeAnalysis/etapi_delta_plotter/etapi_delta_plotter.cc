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

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "AMPTOOLS_DATAIO/EtaPiDeltaPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"

typedef EtaPiDeltaPlotGenerator PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( TwoPSAngles() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
}

//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

int main( int argc, char* argv[] ){
    
    
    // ************************
    // usage
    // ************************
    
    cout << endl << " *** Viewing Results Using AmpPlotter and writing root histograms *** " << endl << endl;
    
    if (argc < 2){
        cout << "Usage:" << endl << endl;
        cout << "\etapi_plotter <results file name> -o <output file name>" << endl << endl;
        return 0;
    }
    
    bool showGui = false;
    string outName = "etapi_plot.root";
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
    //string resultsName(argv[1]);
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
    
    string reactionName = results.reactionList()[0];
    plotGen.enableReaction( reactionName );
    vector<string> sums = plotGen.uniqueSums();
    
    
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
        
        
        // loop over data, accMC, and genMC and kBkgnd
        for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
            if (isum < sums.size() && iplot == PlotGenerator::kData) continue; // only plot data once
            
            // loop over different variables
            for (unsigned int ivar  = 0; ivar  < EtaPiDeltaPlotGenerator::kNumHists; ivar++){
                
                // set unique histogram name for each plot (could put in directories...)
                string histname =  "";
                if (ivar == EtaPiDeltaPlotGenerator::kEtaPiMass)  histname += "EtaPiMass";
                else if (ivar == EtaPiDeltaPlotGenerator::kDeltaPPMass)  histname += "DeltaPPMass";
                else if (ivar == EtaPiDeltaPlotGenerator::kEtaCosTheta)  histname += "cosTheta";
                else if (ivar == EtaPiDeltaPlotGenerator::kPhi)  histname += "Phi";
                else if (ivar == EtaPiDeltaPlotGenerator::kt)  histname += "t";
                else continue;
                
                if (iplot == PlotGenerator::kData) histname += "dat";
                if (iplot == PlotGenerator::kAccMC) histname += "acc";
                if (iplot == PlotGenerator::kGenMC) histname += "gen";
                //if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";
                
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
    // start the GUI
    // ************************
    
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
    
    return 0;
    
}

