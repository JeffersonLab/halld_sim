#include <iostream>
#include <string>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/Flatte.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Ylm.h"
#include "AMPTOOLS_AMPS/Ylm_reim.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_DATAIO/etaetapPlotGeneratorYlm.h"

typedef etaetapPlotGeneratorYlm PlotGen;

void atiSetup(){
  
  // the PlotGenerator will create an AmpToolsInterface in order
  // to create plots - this setup must happen before the
  // AmpToolsInterface is created

  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( Flatte() );
  AmpToolsInterface::registerAmplitude( Ylm());
  AmpToolsInterface::registerAmplitude( Ylm_reim());
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );
}


//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter *** " << endl << endl;

  if (argc < 2 || argc>3){
    
    cout << "Usage:" << endl << endl;
    cout << "\tampPlotter <fit results name> <optional: output pdf name>" << endl << endl;
    return 0;
  }

    // ************************
    // parse the command line parameters
    // ************************

  string resultsName(argv[1]);
  TString outputName = "";
  bool runInteractive=true;
  if(argc==3){
    runInteractive=false;
    outputName=argv[2];
  }

  FitResults results( resultsName );
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }
  string reactionName = results.reactionList()[0];
    // ************************
    // set up the plot generator
    // ************************

  atiSetup();  
  PlotGen plotGen( results );


  if(runInteractive){
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
  }else{
    vector<string> amps = plotGen.uniqueAmplitudes();                                                                 
    // enable all amps for now, eventually want to single out just p wave too                                         
    for(unsigned int i=0;i<amps.size();i++){                                                                          
      plotGen.enableAmp(i);                                                                                           
    }


    // loop over different variables                                                                                  
    TCanvas *c = new TCanvas("c","c",1500,500);                                                                       
    c->Divide(3);                                                                                                     
    gStyle->SetOptStat(0);                                                                                            
    TH1* dataHist;                                                                                                    
    for (unsigned int ivar  = 0; ivar  < PlotGen::kNumHists; ivar++){                                                 
                                                                                                                      
      if(ivar==2||ivar==7||ivar==8){                                                                                  
        cout << "making plot for ivar="<<ivar<<endl;                                                                  
        if(ivar==2) c->cd(1);                                                                                         
        if(ivar==7) c->cd(2);                                                                                         
        if(ivar==8) c->cd(3);                                                                                         
      }else{                                                                                                          
        continue;                                                                                                     
      }                                                                                                               
                                                                                                                      
      for (unsigned int iplot = 0; iplot < PlotGen::kNumTypes; iplot++){                                              
        bool drawFirst = false;                                                                                       
        if(iplot==PlotGen::kData){                                                                                    
          drawFirst = true;                                                                                           
        }else if(iplot!=PlotGen::kAccMC){                                                                             
          continue;                                                                                                   
        }                                                                                                             
                                                                                                                      
        TH1 *hist = plotGen.projection(ivar,reactionName, iplot)->toRoot();                                           
                                                                                                                      
        if(drawFirst){                                                                                                
          hist->Draw();                                                                                               
          dataHist = hist;                                                                                            
        }else{                                                                                                        
          // set color to roughly match GlueX standard                                                                
          hist->SetFillColor(kGreen-8);                                                                               
          hist->Draw("sameHIST");                                                                                     
          dataHist->Draw("same");                                                                                     
          gPad->RedrawAxis();                                                                                         
        }                                                                                                             
                                                                                                                      
      }                                                                                                               
    }                                                                                                                 
    c->Print(outputName);                                                                                             
  }

  cout << "SUCCESS" << endl;
  return 0;

}

