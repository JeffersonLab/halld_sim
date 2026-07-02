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

#include "AMPTOOLS_DATAIO/TwoPiDeltaPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderHist.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta.h"
#include "AMPTOOLS_AMPS/DeltaAngles_Rho.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_only_upperv.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_DoubleSDMEs.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_DoubleSDMEs_unpol.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_factorized.h"
//#include "AMPTOOLS_AMPS/Delta_rho_tdep.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"

typedef TwoPiDeltaPlotGenerator PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_only_upperv() );
  AmpToolsInterface::registerAmplitude( DeltaAngles_Rho() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_DoubleSDMEs() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_DoubleSDMEs_unpol() );
  AmpToolsInterface::registerAmplitude( TwoPiAngles_Delta_factorized() );
  // AmpToolsInterface::registerAmplitude( Delta_rho_tdep() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderHist() );
}

using namespace std;

int main( int argc, char* argv[] ){ //char array for strings


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter and writing root histograms *** " << endl << endl;

  if (argc < 2){
    cout << "Usage:" << endl << endl;
    cout << "\ttwopidelta_plotter <results file name> -o <output file name>" << endl << endl;
    return 0;
  }

  bool showGui = false;
  bool rootfile_creation = true;
  string outName = "twopidelta_plot.root";
  string resultsName(argv[1]);
  Double_t tmean;
  Double_t terr;
  double tmean_old = 0;
  double terr_old = 0;
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
    if (arg == "-n"){
      rootfile_creation = false;
      std::ifstream infile("twopidelta_fitPars.txt");
      std::string name;
      double value, error;
      while (infile >> name >> value >> error) {
          if (name == "tmean") {
              tmean_old = value;
              terr_old = error;
              break;
          }
      }
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

  FitResults results( resultsName ); //get fitresults
  cout << "Fit results loaded" << endl;
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }

    // ************************
    // set up the plot generator
    // ************************
  cout << "Setting up ATI" << endl;
  atiSetup();
  cout << "ATI set up" << endl;
  PlotGen plotGen( results );
  cout << "Plot generator set up" << endl;

    // ************************
    // set up an output ROOT file to store histograms
    // ************************
  if (rootfile_creation){
  TFile* plotfile = new TFile( outName.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  string reactionName = results.reactionList()[0];
  auto reactions = results.reactionList(); 
  cout << "Number of reactions: " << reactions.size() << endl;
  for (size_t k = 0; k < reactions.size(); ++k) {
    auto& reaction = reactions[k];
  string reactionName = reaction;
  cout << "Reaction name:  " << reactionName << endl;
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


    // loop over data, accMC, and genMC
    for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
      if (isum < sums.size() && iplot == PlotGenerator::kData) continue; // only plot data once

      // loop over different variables
      for (unsigned int ivar  = 0; ivar  < TwoPiDeltaPlotGenerator::kNumHists; ivar++){

        // set unique histogram name for each plot (could put in directories...)
        string histname =  "";
        if (ivar == TwoPiDeltaPlotGenerator::k2PiMass) histname += "M2pi";
        else if (ivar == TwoPiDeltaPlotGenerator::kPPipMass) histname += "Mppip";
        else if (ivar == TwoPiDeltaPlotGenerator::kPPimMass) histname += "Mppim";
        else if (ivar == TwoPiDeltaPlotGenerator::kPi0PimMass) histname += "Mpi0pim";
        else if (ivar == TwoPiDeltaPlotGenerator::kPPipPimMass) histname += "Mppippim";
        else if (ivar == TwoPiDeltaPlotGenerator::kPPipPi0Mass) histname += "Mppippi0";
        else if (ivar == TwoPiDeltaPlotGenerator::kPPimPi0Mass) histname += "Mppimpi0";
        else if (ivar == TwoPiDeltaPlotGenerator::kPipPi0PimMass) histname += "Mpippi0pim";
        else if (ivar == TwoPiDeltaPlotGenerator::kPipPimMass) histname += "Mpippim";     
        else if (ivar == TwoPiDeltaPlotGenerator::kPPi0Mass) histname += "Mppi0";           
        else if (ivar == TwoPiDeltaPlotGenerator::kPipPi0Mass) histname += "Mpippi0";      


        else if (ivar == TwoPiDeltaPlotGenerator::kCosTheta_GFJ_Delta) histname += "cosTheta_GFJ_Delta";
        else if (ivar == TwoPiDeltaPlotGenerator::kCosTheta_HF_Delta) histname += "cosTheta_HF_Delta";
        else if (ivar == TwoPiDeltaPlotGenerator::kCosTheta_HF_rho) histname += "cosTheta_HF_rho";
        else if (ivar == TwoPiDeltaPlotGenerator::kCosTheta_GFJ_rho) histname += "cosTheta_GFJ_rho";



        else if (ivar == TwoPiDeltaPlotGenerator::kPhiPiPlus) histname += "PhiPiPlus";
        else if (ivar == TwoPiDeltaPlotGenerator::kPhiPiMinus) histname += "PhiPiMinus";
        else if (ivar == TwoPiDeltaPlotGenerator::kPhiProton) histname += "PhiProton";
        else if (ivar == TwoPiDeltaPlotGenerator::kThetaPiPlus) histname += "ThetaPiPlus";
        else if (ivar == TwoPiDeltaPlotGenerator::kThetaPiMinus) histname += "ThetaPiMinus";
        else if (ivar == TwoPiDeltaPlotGenerator::kThetaDelta) histname += "ThetaDelta";
        else if (ivar == TwoPiDeltaPlotGenerator::kMomPiPlus) histname += "MomPiPlus";
        else if (ivar == TwoPiDeltaPlotGenerator::kMomPiMinus) histname += "MomPiMinus";
        else if (ivar == TwoPiDeltaPlotGenerator::kMomPi0) histname += "MomPi0";
        else if (ivar == TwoPiDeltaPlotGenerator::kMomProton) histname += "MomProton";
        else if (ivar == TwoPiDeltaPlotGenerator::klongMomPiPlus) histname += "longMomPiPlus";
        else if (ivar == TwoPiDeltaPlotGenerator::klongMomPiMinus) histname += "longMomPiMinus";
        else if (ivar == TwoPiDeltaPlotGenerator::klongMomPi0) histname += "longMomPi0";
        else if (ivar == TwoPiDeltaPlotGenerator::klongMomProton) histname += "longMomProton";
        else if (ivar == TwoPiDeltaPlotGenerator::kPhi) histname += "Phi";
        else if (ivar == TwoPiDeltaPlotGenerator::kPhi_GFJ) histname += "PhiGFJ";
        else if (ivar == TwoPiDeltaPlotGenerator::kphi_GFJ_Delta) histname += "phi_GFJ_Delta";
        else if (ivar == TwoPiDeltaPlotGenerator::kphi_HF_Delta) histname += "phi_HF_Delta";
        else if (ivar == TwoPiDeltaPlotGenerator::kphi_HF_rho) histname += "phi_HF_rho";
        else if (ivar == TwoPiDeltaPlotGenerator::kphi_GFJ_rho) histname += "phi_GFJ_rho";
        else if (ivar == TwoPiDeltaPlotGenerator::kPsi) histname += "psi";
        else if (ivar == TwoPiDeltaPlotGenerator::kt) histname += "-t";
        else if (ivar == TwoPiDeltaPlotGenerator::kBeamassymetrie) histname += "Beamass";
        else if (ivar == TwoPiDeltaPlotGenerator::kBeamassymetrie_Delta) histname += "Beamass_Delta";

        else continue;

        if (iplot == PlotGenerator::kData) histname += "dat";
        if (iplot == PlotGenerator::kAccMC) histname += "acc";
        if (iplot == PlotGenerator::kGenMC) histname += "gen";
        if (iplot == PlotGenerator::kBkgnd) histname += "_bkgnd";

        if (isum < sums.size()){
          //ostringstream sdig;  sdig << (isum + 1);
          //histname += sdig.str();

	  // get name of sum for naming histogram
          string sumName = sums[isum];
          histname += "_";
          histname += sumName;
        }
        ;
        if (k==0){
        Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
        TH1* thist = hist->toRoot();
        thist->SetName(histname.c_str());
        
        if (histname=="-tdat"&& k==reactions.size()-1){ 
        tmean = thist->GetMean();
        terr = thist->GetRMS();
          }
        plotfile->cd();
        thist->Write();
        }
        else
        {
        Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
        TH1* histSumReac = (TH1*) plotfile->Get(histname.c_str());
        TH1* thist = hist->toRoot();
        histSumReac->Add(thist);
        if (histname=="-tdat" && k==reactions.size()-1){ 
        tmean = histSumReac->GetMean();
        terr = histSumReac->GetRMS();
        }
        plotfile->cd();
        histSumReac->Write(histname.c_str(), TObject::kOverwrite);
        
        }

        histname += "_";
        histname += reactionName;
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
    // ************************
    // retrieve SDME parameters for plotting and asymmetry
    // ************************

  // parameters to check! be careful! for rho last number is upper case # Delta first number
  vector< string > pars;
  pars.push_back("rho000");
  pars.push_back("rho100");
  pars.push_back("rho1m10");

  pars.push_back("rho111");
  pars.push_back("rho001");
  pars.push_back("rho101");
  pars.push_back("rho1m11");

  pars.push_back("rho102");
  pars.push_back("rho1m12");
	pars.push_back("delta_rho011");
	pars.push_back("delta_rho031");
	pars.push_back("delta_rho03m1");
	pars.push_back("delta_rho111");
	pars.push_back("delta_rho133");
	pars.push_back("delta_rho131");
	pars.push_back("delta_rho13m1");
	pars.push_back("delta_rho231");
	pars.push_back("delta_rho23m1");
  pars.push_back("r00_33_0");
  pars.push_back("r00_11_0");
  pars.push_back("r11_33_0");
  pars.push_back("r11_11_0");
  pars.push_back("r1m1_33_0");
  pars.push_back("r1m1_11_0");
  pars.push_back("r10_33_0");
  pars.push_back("r10_11_0");
  pars.push_back("r00_31_0");
  pars.push_back("r00_3m1_0");
  pars.push_back("r11_31_0");
  pars.push_back("r11_3m1_0");
  pars.push_back("r1m1_31_0");
  pars.push_back("r1m1_3m1_0");
  pars.push_back("r10_31_0");
  pars.push_back("r10_3m1_0");
  pars.push_back("rt1m1_31_0");
  pars.push_back("rt1m1_3m1_0");
  pars.push_back("rt10_31_0");
  pars.push_back("rt10_3m1_0");
  // alpha = 1
  pars.push_back("r00_33_1");
  pars.push_back("r00_11_1");
  pars.push_back("r11_33_1");
  pars.push_back("r11_11_1");
  pars.push_back("r1m1_33_1");
  pars.push_back("r1m1_11_1");
  pars.push_back("r10_33_1");
  pars.push_back("r10_11_1");
  pars.push_back("r00_31_1");
  pars.push_back("r00_3m1_1");
  pars.push_back("r11_31_1");
  pars.push_back("r11_3m1_1");
  pars.push_back("r1m1_31_1");
  pars.push_back("r1m1_3m1_1");
  pars.push_back("r10_31_1");
  pars.push_back("r10_3m1_1");
  pars.push_back("rt1m1_31_1");
  pars.push_back("rt1m1_3m1_1");
  pars.push_back("rt10_31_1");
  pars.push_back("rt10_3m1_1");

  // alpha = 2
  pars.push_back("rt00_31_2");
  pars.push_back("rt00_3m1_2");
  pars.push_back("rt11_31_2");
  pars.push_back("rt11_3m1_2");
  pars.push_back("r1m1_33_2");
  pars.push_back("r1m1_11_2");
  pars.push_back("r1m1_31_2");
  pars.push_back("r1m1_3m1_2");
  pars.push_back("rt1m1_31_2");
  pars.push_back("rt1m1_3m1_2");
  pars.push_back("r10_33_2");
  pars.push_back("r10_11_2");
  pars.push_back("r10_31_2");
  pars.push_back("r10_3m1_2");
  pars.push_back("rt10_31_2");
  pars.push_back("rt10_3m1_2");
  pars.push_back("a");
  pars.push_back("b");
    // pars.push_back("c");
    // pars.push_back("d");
	


  // file for writing parameters (later switch to putting in ROOT file)
  ofstream outfile;
  outfile.open( "twopidelta_fitPars.txt" );
  if (rootfile_creation){
  outfile << "tmean" << "\t" << tmean << "\t" << terr << "\t"  << "\n";
  }
  else{
  outfile << "tmean" << "\t" << tmean_old << "\t" << terr_old << "\t"  << "\n";
  }
  for(unsigned int i = 0; i<pars.size(); i++) {
    
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    if (std::isnan(parValue) || std::isnan(parError)) {
      continue;
  }
  if (pars[i]=="rho001") {
    outfile << "rho001_upperv" << "\t" << parValue << "\t" << parError << "\t"  << "\n";
    continue;  
    }
    outfile << pars[i] << "\t" << parValue << "\t" << parError << "\t"  << "\n";
    
  }
  double rho001_value = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"));
  double rho001_error = 2*sqrt(results.parError("delta_rho111")*results.parError("delta_rho111")+results.parError("delta_rho133")*results.parError("delta_rho133")+results.parError("rho111")*results.parError("rho111")+2*(results.covariance("delta_rho111","delta_rho133")-results.covariance("delta_rho111","rho111")-results.covariance("delta_rho133","rho111")));
  outfile << "rho001" << "\t" <<  rho001_value << "\t" <<  rho001_error << "\t" << "\n";

  double P_sigma_value = 2*results.parValue("rho1m11") - rho001_value; //error calculation uses rho001 expressed by the three other SDMEs
  double P_sigma_error = 2*sqrt(results.parError("rho1m11")*results.parError("rho1m11")+results.parError("rho111")*results.parError("rho111")+results.parError("delta_rho111")*results.parError("delta_rho111")+results.parError("delta_rho133")*results.parError("delta_rho133")+2*(results.covariance("rho1m11","rho111")-results.covariance("rho1m11","delta_rho111")-results.covariance("rho1m11","delta_rho133")-results.covariance("rho111","delta_rho111")-results.covariance("rho111","delta_rho133")+results.covariance("delta_rho111","delta_rho133")));
  outfile << "P" << "\t" <<  P_sigma_value << "\t" <<  P_sigma_error << "\t"  << "\n";




  double r00_33_0_single_value = results.parValue("rho000") * (0.5-results.parValue("delta_rho011"));
  double r00_33_0_single_error = sqrt(results.parError("rho000")*results.parError("rho000")*(0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho000")*results.parValue("rho000")-2*results.covariance("rho000","delta_rho011")*(0.5-results.parValue("delta_rho011"))*results.parValue("rho000"));
  outfile << "r00_33_0_single" << "\t" << r00_33_0_single_value << "\t" << r00_33_0_single_error << "\t" << "\n";

  double r00_31_0_single_value = results.parValue("rho000") * results.parValue("delta_rho031");
  double r00_31_0_single_error = sqrt(results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho031")*results.parValue("delta_rho031")+results.parError("delta_rho031")*results.parError("delta_rho031")*results.parValue("rho000")*results.parValue("rho000")+2*results.covariance("rho000","delta_rho031")*results.parValue("delta_rho031")*results.parValue("rho000"));
  outfile << "r00_31_0_single" << "\t" << r00_31_0_single_value << "\t" << r00_31_0_single_error << "\t" << "\n";

  double r00_11_0_single_value = results.parValue("rho000") * results.parValue("delta_rho011");
  double r00_11_0_single_error = sqrt(results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho011")*results.parValue("delta_rho011")+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho000")*results.parValue("rho000")+2*results.covariance("rho000","delta_rho011")*results.parValue("delta_rho011")*results.parValue("rho000"));
  outfile << "r00_11_0_single" << "\t" << r00_11_0_single_value << "\t" << r00_11_0_single_error << "\t" << "\n";

  double r00_3m1_0_single_value = results.parValue("rho000") * results.parValue("delta_rho03m1");
  double r00_3m1_0_single_error = sqrt(results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")+results.parError("delta_rho03m1")*results.parError("delta_rho03m1")*results.parValue("rho000")*results.parValue("rho000")+2*results.covariance("rho000","delta_rho03m1")*results.parValue("delta_rho03m1")*results.parValue("rho000"));
  outfile << "r00_3m1_0_single" << "\t" << r00_3m1_0_single_value << "\t" << r00_3m1_0_single_error << "\t" << "\n";


  double r11_33_0_single_value = (1-results.parValue("rho000")) * 0.5 * (0.5-results.parValue("delta_rho011"));
  double r11_33_0_single_error = sqrt(0.25*(0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho000")*results.parError("rho000")+0.25*(results.parValue("rho000")-1)*(results.parValue("rho000")-1)*results.parError("delta_rho011")*results.parError("delta_rho011")-0.5*(results.parValue("rho000")-1)*(0.5-results.parValue("delta_rho011"))*results.covariance("rho000","delta_rho011"));
  outfile << "r11_33_0_single" << "\t" << r11_33_0_single_value << "\t" << r11_33_0_single_error << "\t" << "\n";

  double r11_31_0_single_value = (1-results.parValue("rho000")) * 0.5 * results.parValue("delta_rho031");
  double r11_31_0_single_error = sqrt(0.25*results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho031")*results.parValue("delta_rho031")+0.25*results.parError("delta_rho031")*results.parError("delta_rho031")*(results.parValue("rho000")-1)*(results.parValue("rho000")-1)+0.5*results.covariance("rho000","delta_rho031")*results.parValue("delta_rho031")*(results.parValue("rho000")-1));
  outfile << "r11_31_0_single" << "\t" << r11_31_0_single_value << "\t" << r11_31_0_single_error << "\t" << "\n";

  double r11_11_0_single_value = (1-results.parValue("rho000")) * 0.5 * results.parValue("delta_rho011");
  double r11_11_0_single_error = sqrt(0.25*results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho011")*results.parValue("delta_rho011")+0.25*results.parError("delta_rho011")*results.parError("delta_rho011")*(results.parValue("rho000")-1)*(results.parValue("rho000")-1)+0.5*results.covariance("rho000","delta_rho011")*results.parValue("delta_rho011")*(results.parValue("rho000")-1));
  outfile << "r11_11_0_single" << "\t" << r11_11_0_single_value << "\t" << r11_11_0_single_error << "\t" << "\n";

  double r11_3m1_0_single_value = (1-results.parValue("rho000")) * 0.5 * results.parValue("delta_rho03m1");
  double r11_3m1_0_single_error = sqrt(0.25*results.parError("rho000")*results.parError("rho000")*results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")+0.25*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")*(results.parValue("rho000")-1)*(results.parValue("rho000")-1)+0.5*results.covariance("rho000","delta_rho03m1")*results.parValue("delta_rho03m1")*(results.parValue("rho000")-1));
  outfile << "r11_3m1_0_single" << "\t" << r11_3m1_0_single_value << "\t" << r11_3m1_0_single_error << "\t" << "\n";


  double r1m1_33_0_single_value = results.parValue("rho1m10") * (0.5-results.parValue("delta_rho011"));
  double r1m1_33_0_single_error = sqrt(results.parError("rho1m10")*results.parError("rho1m10")*(0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho1m10")*results.parValue("rho1m10")-2*results.covariance("rho1m10","delta_rho011")*(0.5-results.parValue("delta_rho011"))*results.parValue("rho1m10"));
  outfile << "r1m1_33_0_single" << "\t" << r1m1_33_0_single_value << "\t" << r1m1_33_0_single_error << "\t" << "\n";

  double r1m1_11_0_single_value = results.parValue("rho1m10") * results.parValue("delta_rho011");
  double r1m1_11_0_single_error = sqrt(results.parError("rho1m10")*results.parError("rho1m10")*results.parValue("delta_rho011")*results.parValue("delta_rho011")+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho1m10")*results.parValue("rho1m10")+2*results.covariance("rho1m10","delta_rho011")*results.parValue("delta_rho011")*results.parValue("rho1m10"));
  outfile << "r1m1_11_0_single" << "\t" << r1m1_11_0_single_value << "\t" << r1m1_11_0_single_error << "\t" << "\n";

  double r1m1_31_0_single_value = results.parValue("rho1m10") * results.parValue("delta_rho031");
  double r1m1_31_0_single_error = sqrt(results.parError("rho1m10")*results.parError("rho1m10")*results.parValue("delta_rho031")*results.parValue("delta_rho031")+results.parError("delta_rho031")*results.parError("delta_rho031")*results.parValue("rho1m10")*results.parValue("rho1m10")+2*results.covariance("rho1m10","delta_rho031")*results.parValue("delta_rho031")*results.parValue("rho1m10"));
  outfile << "r1m1_31_0_single" << "\t" << r1m1_31_0_single_value << "\t" << r1m1_31_0_single_error << "\t" << "\n";

  double r1m1_3m1_0_single_value = results.parValue("rho1m10") * results.parValue("delta_rho03m1");
  double r1m1_3m1_0_single_error = sqrt(results.parError("rho1m10")*results.parError("rho1m10")*results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")+results.parError("delta_rho03m1")*results.parError("delta_rho03m1")*results.parValue("rho1m10")*results.parValue("rho1m10")+2*results.covariance("rho1m10","delta_rho03m1")*results.parValue("delta_rho03m1")*results.parValue("rho1m10"));
  outfile << "r1m1_3m1_0_single" << "\t" << r1m1_3m1_0_single_value << "\t" << r1m1_3m1_0_single_error << "\t" << "\n";


  double r10_33_0_single_value = results.parValue("rho100") * (0.5-results.parValue("delta_rho011"));
  double r10_33_0_single_error = sqrt(results.parError("rho100")*results.parError("rho100")*(0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho100")*results.parValue("rho100")-2*results.covariance("rho100","delta_rho011")*(0.5-results.parValue("delta_rho011"))*results.parValue("rho100"));
  outfile << "r10_33_0_single" << "\t" << r10_33_0_single_value << "\t" << r10_33_0_single_error << "\t" << "\n";

  double r10_31_0_single_value = results.parValue("rho100") * results.parValue("delta_rho031");
  double r10_31_0_single_error = sqrt(results.parError("rho100")*results.parError("rho100")*results.parValue("delta_rho031")*results.parValue("delta_rho031")+results.parError("delta_rho031")*results.parError("delta_rho031")*results.parValue("rho100")*results.parValue("rho100")+2*results.covariance("rho100","delta_rho031")*results.parValue("delta_rho031")*results.parValue("rho100"));
  outfile << "r10_31_0_single" << "\t" << r10_31_0_single_value << "\t" << r10_31_0_single_error << "\t" << "\n";

  double r10_11_0_single_value = results.parValue("rho100") * results.parValue("delta_rho011");
  double r10_11_0_single_error = sqrt(results.parError("rho100")*results.parError("rho100")*results.parValue("delta_rho011")*results.parValue("delta_rho011")+results.parError("delta_rho011")*results.parError("delta_rho011")*results.parValue("rho100")*results.parValue("rho100")+2*results.covariance("rho100","delta_rho011")*results.parValue("delta_rho011")*results.parValue("rho100"));
  outfile << "r10_11_0_single" << "\t" << r10_11_0_single_value << "\t" << r10_11_0_single_error << "\t" << "\n";

  double r10_3m1_0_single_value = results.parValue("rho100") * results.parValue("delta_rho03m1");
  double r10_3m1_0_single_error = sqrt(results.parError("rho100")*results.parError("rho100")*results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")+results.parError("delta_rho03m1")*results.parError("delta_rho03m1")*results.parValue("rho100")*results.parValue("rho100")+2*results.covariance("rho100","delta_rho03m1")*results.parValue("delta_rho03m1")*results.parValue("rho100"));
  outfile << "r10_3m1_0_single" << "\t" << r10_3m1_0_single_value << "\t" << r10_3m1_0_single_error << "\t" << "\n";



  //alpha=1
  double r00_33_1_single_value = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")) * results.parValue("delta_rho133");
  double r00_33_1_single_error = sqrt((results.parValue("delta_rho133")*results.parValue("delta_rho133"))*4*((results.parError("delta_rho111")*results.parError("delta_rho111"))+(results.parError("delta_rho133")*results.parError("delta_rho133"))+(results.parError("rho111")*results.parError("rho111"))+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(results.parError("delta_rho133")*results.parError("delta_rho133"))+2*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*results.parValue("delta_rho133")*2*(results.covariance("delta_rho111","delta_rho133")+results.covariance("delta_rho133","delta_rho133")-results.covariance("rho111","delta_rho133")));  
  outfile << "r00_33_1_single" << "\t" << r00_33_1_single_value << "\t" << r00_33_1_single_error << "\n";
  
  double r00_31_1_single_value = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")) * results.parValue("delta_rho131");
  double r00_31_1_single_error = sqrt((results.parValue("delta_rho131")*results.parValue("delta_rho131"))*4*((results.parError("delta_rho111")*results.parError("delta_rho111"))+(results.parError("delta_rho133")*results.parError("delta_rho133"))+(results.parError("rho111")*results.parError("rho111"))+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(results.parError("delta_rho131")*results.parError("delta_rho131"))+2*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*results.parValue("delta_rho131")*2*(results.covariance("delta_rho111","delta_rho131")+results.covariance("delta_rho133","delta_rho131")-results.covariance("rho111","delta_rho131")));
  outfile << "r00_31_1_single" << "\t" << r00_31_1_single_value << "\t" << r00_31_1_single_error << "\t" << "\n";

  double r00_11_1_single_value = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")) * results.parValue("delta_rho111");
  double r00_11_1_single_error = sqrt((results.parValue("delta_rho111")*results.parValue("delta_rho111"))*4*((results.parError("delta_rho111")*results.parError("delta_rho111"))+(results.parError("delta_rho133")*results.parError("delta_rho133"))+(results.parError("rho111")*results.parError("rho111"))+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(results.parError("delta_rho111")*results.parError("delta_rho111"))+2*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*results.parValue("delta_rho111")*2*(results.covariance("delta_rho111","delta_rho111")+results.covariance("delta_rho133","delta_rho111")-results.covariance("rho111","delta_rho111")));
  outfile << "r00_11_1_single" << "\t" << r00_11_1_single_value << "\t" << r00_11_1_single_error << "\t" << "\n";

  double r00_3m1_1_single_value = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")) * results.parValue("delta_rho13m1");
  double r00_3m1_1_single_error = sqrt((results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1"))*4*((results.parError("delta_rho111")*results.parError("delta_rho111"))+(results.parError("delta_rho133")*results.parError("delta_rho133"))+(results.parError("rho111")*results.parError("rho111"))+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*(results.parError("delta_rho13m1")*results.parError("delta_rho13m1"))+2*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111")))*results.parValue("delta_rho13m1")*2*(results.covariance("delta_rho111","delta_rho13m1")+results.covariance("delta_rho133","delta_rho13m1")-results.covariance("rho111","delta_rho13m1")));
  outfile << "r00_3m1_1_single" << "\t" << r00_3m1_1_single_value << "\t" << r00_3m1_1_single_error << "\n";
  
  double r11_33_1_single_value = results.parValue("rho111") * results.parValue("delta_rho133");
  double r11_33_1_single_error = sqrt(results.parError("rho111")*results.parError("rho111")*results.parValue("delta_rho133")*results.parValue("delta_rho133")+results.parError("delta_rho133")*results.parError("delta_rho133")*results.parValue("rho111")*results.parValue("rho111")+2*results.covariance("rho111","delta_rho133")*results.parValue("delta_rho133")*results.parValue("rho111"));
  outfile << "r11_33_1_single" << "\t" << r11_33_1_single_value << "\t" << r11_33_1_single_error << "\n";
  
  double r11_31_1_single_value = results.parValue("rho111") * results.parValue("delta_rho131");
  double r11_31_1_single_error = sqrt(results.parError("rho111")*results.parError("rho111")*results.parValue("delta_rho131")*results.parValue("delta_rho131")+results.parError("delta_rho131")*results.parError("delta_rho131")*results.parValue("rho111")*results.parValue("rho111")+2*results.covariance("rho111","delta_rho131")*results.parValue("delta_rho131")*results.parValue("rho111"));
  outfile << "r11_31_1_single" << "\t" << r11_31_1_single_value << "\t" << r11_31_1_single_error << "\n";
  
  double r11_11_1_single_value = results.parValue("rho111") * results.parValue("delta_rho111");
  double r11_11_1_single_error = sqrt(results.parError("rho111")*results.parError("rho111")*results.parValue("delta_rho111")*results.parValue("delta_rho111")+results.parError("delta_rho111")*results.parError("delta_rho111")*results.parValue("rho111")*results.parValue("rho111")+2*results.covariance("rho111","delta_rho111")*results.parValue("delta_rho111")*results.parValue("rho111"));
  outfile << "r11_11_1_single" << "\t" << r11_11_1_single_value << "\t" << r11_11_1_single_error << "\n";
  
  double r11_3m1_1_single_value = results.parValue("rho111") * results.parValue("delta_rho13m1");
  double r11_3m1_1_single_error = sqrt(results.parError("rho111")*results.parError("rho111")*results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")+results.parError("delta_rho13m1")*results.parError("delta_rho13m1")*results.parValue("rho111")*results.parValue("rho111")+2*results.covariance("rho111","delta_rho13m1")*results.parValue("delta_rho13m1")*results.parValue("rho111"));
  outfile << "r11_3m1_1_single" << "\t" << r11_3m1_1_single_value << "\t" << r11_3m1_1_single_error << "\n";
  
  double r1m1_33_1_single_value = results.parValue("rho1m11") * results.parValue("delta_rho133");
  double r1m1_33_1_single_error = sqrt(results.parError("rho1m11")*results.parError("rho1m11")*results.parValue("delta_rho133")*results.parValue("delta_rho133")+results.parError("delta_rho133")*results.parError("delta_rho133")*results.parValue("rho1m11")*results.parValue("rho1m11")+2*results.covariance("rho1m11","delta_rho133")*results.parValue("delta_rho133")*results.parValue("rho1m11"));
  outfile << "r1m1_33_1_single" << "\t" << r1m1_33_1_single_value << "\t" << r1m1_33_1_single_error << "\n";
  
  double r1m1_11_1_single_value = results.parValue("rho1m11") * results.parValue("delta_rho111");
  double r1m1_11_1_single_error = sqrt(results.parError("rho1m11")*results.parError("rho1m11")*results.parValue("delta_rho111")*results.parValue("delta_rho111")+results.parError("delta_rho111")*results.parError("delta_rho111")*results.parValue("rho1m11")*results.parValue("rho1m11")+2*results.covariance("rho1m11","delta_rho111")*results.parValue("delta_rho111")*results.parValue("rho1m11"));
  outfile << "r1m1_11_1_single" << "\t" << r1m1_11_1_single_value << "\t" << r1m1_11_1_single_error << "\n";
  
  double r1m1_31_1_single_value = results.parValue("rho1m11") * results.parValue("delta_rho131");
  double r1m1_31_1_single_error = sqrt(results.parError("rho1m11")*results.parError("rho1m11")*results.parValue("delta_rho131")*results.parValue("delta_rho131")+results.parError("delta_rho131")*results.parError("delta_rho131")*results.parValue("rho1m11")*results.parValue("rho1m11")+2*results.covariance("rho1m11","delta_rho131")*results.parValue("delta_rho131")*results.parValue("rho1m11"));
  outfile << "r1m1_31_1_single" << "\t" << r1m1_31_1_single_value << "\t" << r1m1_31_1_single_error << "\n";
  
  double r1m1_3m1_1_single_value = results.parValue("rho1m11") * results.parValue("delta_rho13m1");
  double r1m1_3m1_1_single_error = sqrt(results.parError("rho1m11")*results.parError("rho1m11")*results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")+results.parError("delta_rho13m1")*results.parError("delta_rho13m1")*results.parValue("rho1m11")*results.parValue("rho1m11")+2*results.covariance("rho1m11","delta_rho13m1")*results.parValue("delta_rho13m1")*results.parValue("rho1m11"));
  outfile << "r1m1_3m1_1_single" << "\t" << r1m1_3m1_1_single_value << "\t" << r1m1_3m1_1_single_error << "\n";
  
  double r10_33_1_single_value = results.parValue("rho101") * results.parValue("delta_rho133");
  double r10_33_1_single_error = sqrt(results.parError("rho101")*results.parError("rho101")*results.parValue("delta_rho133")*results.parValue("delta_rho133")+results.parError("delta_rho133")*results.parError("delta_rho133")*results.parValue("rho101")*results.parValue("rho101")+2*results.covariance("rho101","delta_rho133")*results.parValue("delta_rho133")*results.parValue("rho101"));
  outfile << "r10_33_1_single" << "\t" << r10_33_1_single_value << "\t" << r10_33_1_single_error << "\n";
  
  double r10_31_1_single_value = results.parValue("rho101") * results.parValue("delta_rho131");
  double r10_31_1_single_error = sqrt(results.parError("rho101")*results.parError("rho101")*results.parValue("delta_rho131")*results.parValue("delta_rho131")+results.parError("delta_rho131")*results.parError("delta_rho131")*results.parValue("rho101")*results.parValue("rho101")+2*results.covariance("rho101","delta_rho131")*results.parValue("delta_rho131")*results.parValue("rho101"));
  outfile << "r10_31_1_single" << "\t" << r10_31_1_single_value << "\t" << r10_31_1_single_error << "\n";
  
  double r10_11_1_single_value = results.parValue("rho101") * results.parValue("delta_rho111");
  double r10_11_1_single_error = sqrt(results.parError("rho101")*results.parError("rho101")*results.parValue("delta_rho111")*results.parValue("delta_rho111")+results.parError("delta_rho111")*results.parError("delta_rho111")*results.parValue("rho101")*results.parValue("rho101")+2*results.covariance("rho101","delta_rho111")*results.parValue("delta_rho111")*results.parValue("rho101"));
  outfile << "r10_11_1_single" << "\t" << r10_11_1_single_value << "\t" << r10_11_1_single_error << "\n";
  
  double r10_3m1_1_single_value = results.parValue("rho101") * results.parValue("delta_rho13m1");
  double r10_3m1_1_single_error = sqrt(results.parError("rho101")*results.parError("rho101")*results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")+results.parError("delta_rho13m1")*results.parError("delta_rho13m1")*results.parValue("rho101")*results.parValue("rho101")+2*results.covariance("rho101","delta_rho13m1")*results.parValue("delta_rho13m1")*results.parValue("rho101"));
  outfile << "r10_3m1_1_single" << "\t" << r10_3m1_1_single_value << "\t" << r10_3m1_1_single_error << "\n";


  //alpha=0

  // TraceB(rho_mm'_lambdalambda')
  double rho000_double_value = 1.0-4.0*results.parValue("r11_11_0") - 4.0 *results.parValue("r11_33_0");
  double rho000_double_error = sqrt(16.0*results.parError("r11_11_0")*results.parError("r11_11_0")+16.0*results.parError("r11_33_0")*results.parError("r11_33_0")+32.0*results.covariance("r11_11_0","r11_33_0"));
  outfile << "rho000_double" << "\t" << rho000_double_value << "\t" << rho000_double_error << "\t" << "\n";

  double rho100_double_value = 2.0*results.parValue("r10_11_0") +2.0*results.parValue("r10_33_0");
  double rho100_double_error = sqrt(4.0*results.parError("r10_11_0")*results.parError("r10_11_0")+4.0*results.parError("r10_33_0")*results.parError("r10_33_0")+8.0*results.covariance("r10_11_0","r10_33_0"));
  outfile << "rho100_double" << "\t" << rho100_double_value << "\t" << rho100_double_error << "\t" << "\n";

  double rho1m10_double_value = 2.0*results.parValue("r1m1_11_0") +2.0*results.parValue("r1m1_33_0");
  double rho1m10_double_error = sqrt(4.0*results.parError("r1m1_11_0")*results.parError("r1m1_11_0")+4.0*results.parError("r1m1_33_0")*results.parError("r1m1_33_0")+8.0*results.covariance("r1m1_11_0","r1m1_33_0"));
  outfile << "rho1m10_double" << "\t" << rho1m10_double_value << "\t" << rho1m10_double_error << "\t" << "\n";

  // TraceM(rho_mm'_lambdalambda')
  double delta_rho011_double_value = 0.5-results.parValue("r00_33_0")-2*results.parValue("r11_33_0");
  double delta_rho011_double_error = sqrt(results.parError("r00_33_0")*results.parError("r00_33_0")+4*results.parError("r11_33_0")*results.parError("r11_33_0")+4*results.covariance("r00_33_0","r11_33_0"));
  outfile << "delta_rho011_double" << "\t" << delta_rho011_double_value << "\t" << delta_rho011_double_error << "\t" << "\n";

  double delta_rho031_double_value = results.parValue("r00_31_0") + 2 * results.parValue("r11_31_0");
  double delta_rho031_double_error = sqrt(results.parError("r00_31_0")*results.parError("r00_31_0")+4*results.parError("r11_31_0")*results.parError("r11_31_0")+4*results.covariance("r00_31_0","r11_31_0"));
  outfile << "delta_rho031_double" << "\t" << delta_rho031_double_value << "\t" << delta_rho031_double_error << "\t" << "\n";

  double delta_rho03m1_double_value = results.parValue("r00_3m1_0") + 2 * results.parValue("r11_3m1_0");
  double delta_rho03m1_double_error = sqrt(results.parError("r00_3m1_0")*results.parError("r00_3m1_0")+4*results.parError("r11_3m1_0")*results.parError("r11_3m1_0")+4*results.covariance("r00_3m1_0","r11_3m1_0"));
  outfile << "delta_rho03m1_double" << "\t" << delta_rho03m1_double_value << "\t" << delta_rho03m1_double_error << "\t" << "\n";

  // alpha = 1
  double rho111_double_value = 2*(results.parValue("r11_33_1") + results.parValue("r11_11_1"));
  double rho111_double_error = sqrt(4*results.parError("r11_33_1")*results.parError("r11_33_1")+4*results.parError("r11_11_1")*results.parError("r11_11_1")+8*results.covariance("r11_33_1","r11_11_1"));
  outfile << "rho111_double" << "\t" << rho111_double_value << "\t" << rho111_double_error << "\t" << "\n";

  double rho1m11_double_value = 2*(results.parValue("r1m1_33_1") + results.parValue("r1m1_11_1"));
  double rho1m11_double_error = sqrt(4*results.parError("r1m1_33_1")*results.parError("r1m1_33_1")+4*results.parError("r1m1_11_1")*results.parError("r1m1_11_1")+8*results.covariance("r1m1_33_1","r1m1_11_1"));
  outfile << "rho1m11_double" << "\t" << rho1m11_double_value << "\t" << rho1m11_double_error << "\t" << "\n";

  double rho101_double_value = 2*(results.parValue("r10_33_1") + results.parValue("r10_11_1"));
  double rho101_double_error = sqrt(4*results.parError("r10_33_1")*results.parError("r10_33_1")+4*results.parError("r10_11_1")*results.parError("r10_11_1")+8*results.covariance("r10_33_1","r10_11_1"));
  outfile << "rho101_double" << "\t" << rho101_double_value << "\t" << rho101_double_error << "\t" << "\n";

  double delta_rho111_double_value = results.parValue("r00_11_1")+2*results.parValue("r11_11_1");
  double delta_rho111_double_error = sqrt(results.parError("r00_11_1")*results.parError("r00_11_1")+4*results.parError("r11_11_1")*results.parError("r11_11_1")+4*results.covariance("r00_11_1","r11_11_1"));
  outfile << "delta_rho111_double" << "\t" << delta_rho111_double_value << "\t" << delta_rho111_double_error << "\t" << "\n";

  double delta_rho133_double_value = results.parValue("r00_33_1")+2*results.parValue("r11_33_1");
  double delta_rho133_double_error = sqrt(results.parError("r00_33_1")*results.parError("r00_33_1")+4*results.parError("r11_33_1")*results.parError("r11_33_1")+4*results.covariance("r00_33_1","r11_33_1"));
  outfile << "delta_rho133_double" << "\t" << delta_rho133_double_value << "\t" << delta_rho133_double_error << "\t" << "\n";

  double delta_rho131_double_value = results.parValue("r00_31_1")+2*results.parValue("r11_31_1");
  double delta_rho131_double_error = sqrt(results.parError("r00_31_1")*results.parError("r00_31_1")+4*results.parError("r11_3m1_1")*results.parError("r11_3m1_1")+4*results.covariance("r00_31_1","r11_3m1_1"));
  outfile << "delta_rho131_double" << "\t" << delta_rho131_double_value << "\t" << delta_rho131_double_error << "\t" << "\n";

  double delta_rho13m1_double_value = results.parValue("r00_3m1_1")+2*results.parValue("r11_3m1_1");
  double delta_rho13m1_double_error = sqrt(results.parError("r00_3m1_1")*results.parError("r00_3m1_1")+4*results.parError("r11_31_1")*results.parError("r11_31_1")+4*results.covariance("r00_3m1_1","r11_31_1"));
  outfile << "delta_rho13m1_double" << "\t" << delta_rho13m1_double_value << "\t" << delta_rho13m1_double_error << "\t" << "\n";
  
  double rho001_double_value = 2.0*(results.parValue("r00_33_1") + results.parValue("r00_11_1"));
  double rho001_double_error = sqrt(4.0*(results.parError("r00_33_1")*results.parError("r00_33_1") + results.parError("r00_11_1")*results.parError("r00_11_1") + 2.0*results.covariance("r00_33_1","r00_11_1")));
  outfile << "rho001_double" << "\t" << rho001_double_value << "\t" << rho001_double_error << "\n";
  // alpha = 2
  double rho102_double_value = 2*(results.parValue("r10_33_2") + results.parValue("r10_11_2"));
  double rho102_double_error = sqrt(4*results.parError("r10_33_2")*results.parError("r10_33_2")+4*results.parError("r10_11_2")*results.parError("r10_11_2")+8*results.covariance("r10_33_2","r10_11_2"));
  outfile << "rho102_double" << "\t" << rho102_double_value << "\t" << rho102_double_error << "\t" << "\n";

  double rho1m12_double_value = 2*(results.parValue("r1m1_33_2") + results.parValue("r1m1_11_2"));
  double rho1m12_double_error = sqrt(4*results.parError("r1m1_33_2")*results.parError("r1m1_33_2")+4*results.parError("r1m1_11_2")*results.parError("r1m1_11_2")+8*results.covariance("r1m1_33_2","r1m1_11_2"));
  outfile << "rho1m12_double" << "\t" << rho1m12_double_value << "\t" << rho1m12_double_error << "\t" << "\n";

  double delta_rho231_double_value = results.parValue("rt00_31_2")+2*results.parValue("rt11_3m1_2");
  double delta_rho231_double_error = sqrt(results.parError("rt00_31_2")*results.parError("rt00_31_2")+4*results.parError("rt11_3m1_2")*results.parError("rt11_3m1_2")+4*results.covariance("rt00_31_2","rt11_3m1_2"));
  outfile << "delta_rho231_double" << "\t" << delta_rho231_double_value << "\t" << delta_rho231_double_error << "\t" << "\n";

  double delta_rho23m1_double_value = results.parValue("rt00_3m1_2")+2*results.parValue("rt11_3m1_2");
  double delta_rho23m1_double_error = sqrt(results.parError("rt00_3m1_2")*results.parError("rt00_3m1_2")+4*results.parError("rt11_3m1_2")*results.parError("rt11_3m1_2")+4*results.covariance("rt00_3m1_2","rt11_3m1_2"));
  outfile << "delta_rho23m1_double" << "\t" << delta_rho23m1_double_value << "\t" << delta_rho23m1_double_error << "\t" << "\n";


  //test
  //test
  double r00_33_1_single_value_test = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*(0.5-results.parValue("delta_rho011")) + results.parValue("rho000")*results.parValue("delta_rho133");
  double r00_33_1_single_error_test = sqrt((2*(0.5-results.parValue("delta_rho011"))*results.parError("delta_rho111"))*(2*(0.5-results.parValue("delta_rho011"))*results.parError("delta_rho111")) + ((2*(0.5-results.parValue("delta_rho011"))+results.parValue("rho000"))*results.parError("delta_rho133"))*((2*(0.5-results.parValue("delta_rho011"))+results.parValue("rho000"))*results.parError("delta_rho133")) + (-2*(0.5-results.parValue("delta_rho011"))*results.parError("rho111"))*(-2*(0.5-results.parValue("delta_rho011"))*results.parError("rho111")) + (-2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho011"))*(-2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho011")) + (results.parValue("delta_rho133")*results.parError("rho000"))*(results.parValue("delta_rho133")*results.parError("rho000")));
  outfile << "r00_33_1_single_test" << "\t" << r00_33_1_single_value_test << "\t" << r00_33_1_single_error_test << "\n";
  
  double r00_31_1_single_value_test = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parValue("delta_rho031") + results.parValue("rho000")*results.parValue("delta_rho131");
  double r00_31_1_single_error_test = sqrt((2*results.parValue("delta_rho031")*results.parError("delta_rho111"))*(2*results.parValue("delta_rho031")*results.parError("delta_rho111")) + (2*results.parValue("delta_rho031")*results.parError("delta_rho133"))*(2*results.parValue("delta_rho031")*results.parError("delta_rho133")) + (-2*results.parValue("delta_rho031")*results.parError("rho111"))*(-2*results.parValue("delta_rho031")*results.parError("rho111")) + (2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho031"))*(2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho031")) + (results.parValue("delta_rho131")*results.parError("rho000"))*(results.parValue("delta_rho131")*results.parError("rho000")) + (results.parValue("rho000")*results.parError("delta_rho131"))*(results.parValue("rho000")*results.parError("delta_rho131")));
  outfile << "r00_31_1_single_test" << "\t" << r00_31_1_single_value_test << "\t" << r00_31_1_single_error_test << "\n";


  double r00_11_1_single_value_test = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parValue("delta_rho011")+results.parValue("rho000")*results.parValue("delta_rho111");
  double r00_11_1_single_error_test = sqrt(4*results.parValue("delta_rho011")*results.parValue("delta_rho011")*(results.parError("delta_rho111")*results.parError("delta_rho111")+results.parError("delta_rho133")*results.parError("delta_rho133")+results.parError("rho111")*results.parError("rho111")+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+4*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho011")*results.parError("delta_rho011")+results.parValue("delta_rho111")*results.parValue("delta_rho111")*results.parError("rho000")*results.parError("rho000")+results.parValue("rho000")*results.parValue("rho000")*results.parError("delta_rho111")*results.parError("delta_rho111")+2*results.parValue("rho000")*results.parValue("delta_rho111")*results.covariance("rho000","delta_rho111"));
  outfile << "r00_11_1_single_test" << "\t" << r00_11_1_single_value_test << "\t" << r00_11_1_single_error_test << "\n";
  
  double r00_3m1_1_single_value_test = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parValue("delta_rho03m1") + results.parValue("rho000")*results.parValue("delta_rho13m1");
  double r00_3m1_1_single_error_test = sqrt(4*results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*(results.parError("delta_rho111")*results.parError("delta_rho111")+results.parError("delta_rho133")*results.parError("delta_rho133")+results.parError("rho111")*results.parError("rho111")+2*results.covariance("delta_rho111","delta_rho133")-2*results.covariance("delta_rho111","rho111")-2*results.covariance("delta_rho133","rho111"))+4*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*(results.parValue("delta_rho111")+results.parValue("delta_rho133")-results.parValue("rho111"))*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")*results.parError("rho000")*results.parError("rho000")+results.parValue("rho000")*results.parValue("rho000")*results.parError("delta_rho13m1")*results.parError("delta_rho13m1")+2*results.parValue("rho000")*results.parValue("delta_rho13m1")*results.covariance("rho000","delta_rho13m1"));
  outfile << "r00_3m1_1_single_test" << "\t" << r00_3m1_1_single_value_test << "\t" << r00_3m1_1_single_error_test << "\n";
  
  double r11_33_1_single_value_test = results.parValue("rho111")*(0.5-results.parValue("delta_rho011"))+(1-results.parValue("rho000"))*0.5*results.parValue("delta_rho133");
  double r11_33_1_single_error_test = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho011")*results.parError("delta_rho011")+0.25*results.parValue("delta_rho133")*results.parValue("delta_rho133")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho133")*results.parError("delta_rho133")-2*results.parValue("rho111")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho111","delta_rho011")-(0.5-results.parValue("delta_rho011"))*results.parValue("delta_rho133")*results.covariance("rho111","rho000")+results.parValue("rho111")*(1-results.parValue("rho000"))*results.covariance("rho111","delta_rho133")+results.parValue("rho111")*results.parValue("delta_rho133")*results.covariance("delta_rho011","rho000")-results.parValue("rho111")*(1-results.parValue("rho000"))*results.covariance("delta_rho011","delta_rho133")-0.5*(1-results.parValue("rho000"))*results.parValue("delta_rho133")*results.covariance("rho000","delta_rho133"));
  outfile << "r11_33_1_single_test" << "\t" << r11_33_1_single_value_test << "\t" << r11_33_1_single_error_test << "\n";
  
  double r11_31_1_single_value_test = results.parValue("rho111")*results.parValue("delta_rho031")+(1-results.parValue("rho000"))*0.5*results.parValue("delta_rho131");
  double r11_31_1_single_error_test = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho031")*results.parError("delta_rho031")+0.25*results.parValue("delta_rho131")*results.parValue("delta_rho131")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho131")*results.parError("delta_rho131")+2*results.parValue("delta_rho031")*results.parValue("rho111")*results.covariance("rho111","delta_rho031")-results.parValue("delta_rho031")*results.parValue("delta_rho131")*results.covariance("rho111","rho000")+results.parValue("delta_rho031")*(1-results.parValue("rho000"))*results.covariance("rho111","delta_rho131")-results.parValue("rho111")*results.parValue("delta_rho131")*results.covariance("delta_rho031","rho000")+results.parValue("rho111")*(1-results.parValue("rho000"))*results.covariance("delta_rho031","delta_rho131")-(1-results.parValue("rho000"))*results.parValue("delta_rho131")*results.covariance("rho000","delta_rho131"));
  outfile << "r11_31_1_single_test" << "\t" << r11_31_1_single_value_test << "\t" << r11_31_1_single_error_test << "\n";
  
  double r11_11_1_single_value_test = results.parValue("rho111")*results.parValue("delta_rho011")+(1-results.parValue("rho000"))*0.5*results.parValue("delta_rho111");
  double r11_11_1_single_error_test = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho011")*results.parError("delta_rho011")+0.25*results.parValue("delta_rho111")*results.parValue("delta_rho111")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho111")*results.parError("delta_rho111")+2*results.parValue("delta_rho011")*results.parValue("rho111")*results.covariance("rho111","delta_rho011")-results.parValue("delta_rho011")*results.parValue("delta_rho111")*results.covariance("rho111","rho000")+results.parValue("delta_rho011")*(1-results.parValue("rho000"))*results.covariance("rho111","delta_rho111")-results.parValue("rho111")*results.parValue("delta_rho111")*results.covariance("delta_rho011","rho000")+results.parValue("rho111")*(1-results.parValue("rho000"))*results.covariance("delta_rho011","delta_rho111")-(1-results.parValue("rho000"))*results.parValue("delta_rho111")*results.covariance("rho000","delta_rho111"));
  outfile << "r11_11_1_single_test" << "\t" << r11_11_1_single_value_test << "\t" << r11_11_1_single_error_test << "\n";
  
  double r11_3m1_1_single_value_test = results.parValue("rho111")*results.parValue("delta_rho03m1")+(1-results.parValue("rho000"))*0.5*results.parValue("delta_rho13m1");
  double r11_3m1_1_single_error_test = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+0.25*results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho13m1")*results.parError("delta_rho13m1")+2*results.parValue("delta_rho03m1")*results.parValue("rho111")*results.covariance("rho111","delta_rho03m1")-results.parValue("delta_rho03m1")*results.parValue("delta_rho13m1")*results.covariance("rho111","rho000")+results.parValue("delta_rho03m1")*(1-results.parValue("rho000"))*results.covariance("rho111","delta_rho13m1")-results.parValue("rho111")*results.parValue("delta_rho13m1")*results.covariance("delta_rho03m1","rho000")+results.parValue("rho111")*(1-results.parValue("rho000"))*results.covariance("delta_rho03m1","delta_rho13m1")-(1-results.parValue("rho000"))*results.parValue("delta_rho13m1")*results.covariance("rho000","delta_rho13m1"));
  outfile << "r11_3m1_1_single_test" << "\t" << r11_3m1_1_single_value_test << "\t" << r11_3m1_1_single_error_test << "\n";
  
  double r1m1_33_1_single_value_test = results.parValue("rho1m11")*(0.5-results.parValue("delta_rho011"))+results.parValue("rho1m10")*results.parValue("delta_rho133");
  double r1m1_33_1_single_error_test = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho011")*results.parError("delta_rho011")+results.parValue("delta_rho133")*results.parValue("delta_rho133")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho133")*results.parError("delta_rho133")-2*results.parValue("rho1m11")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho1m11","delta_rho011")+2*(0.5-results.parValue("delta_rho011"))*results.parValue("delta_rho133")*results.covariance("rho1m11","rho1m10")+2*results.parValue("rho1m11")*results.parValue("rho1m10")*results.covariance("rho1m11","delta_rho133")-2*results.parValue("rho1m11")*results.parValue("delta_rho133")*results.covariance("delta_rho011","rho1m10")-2*results.parValue("rho1m11")*results.parValue("rho1m10")*results.covariance("delta_rho011","delta_rho133")+2*results.parValue("rho1m10")*results.parValue("delta_rho133")*results.covariance("rho1m10","delta_rho133"));
  outfile << "r1m1_33_1_single_test" << "\t" << r1m1_33_1_single_value_test << "\t" << r1m1_33_1_single_error_test << "\n";
  
  double r1m1_11_1_single_value_test = results.parValue("rho1m11")*results.parValue("delta_rho011")+results.parValue("rho1m10")*results.parValue("delta_rho111");
  double r1m1_11_1_single_error_test = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho011")*results.parError("delta_rho011")+results.parValue("delta_rho111")*results.parValue("delta_rho111")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho111")*results.parError("delta_rho111")+2*results.parValue("delta_rho011")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho011")+2*results.parValue("delta_rho011")*results.parValue("delta_rho111")*results.covariance("rho1m11","rho1m10")+2*results.parValue("delta_rho011")*results.parValue("rho1m10")*results.covariance("rho1m11","delta_rho111")+2*results.parValue("rho1m11")*results.parValue("delta_rho111")*results.covariance("delta_rho011","rho1m10")+2*results.parValue("rho1m11")*results.parValue("rho1m10")*results.covariance("delta_rho011","delta_rho111")+2*results.parValue("rho1m10")*results.parValue("delta_rho111")*results.covariance("rho1m10","delta_rho111"));
  outfile << "r1m1_11_1_single_test" << "\t" << r1m1_11_1_single_value_test << "\t" << r1m1_11_1_single_error_test << "\n";
  
  double r1m1_31_1_single_value_test = results.parValue("rho1m11")*results.parValue("delta_rho031")+results.parValue("rho1m10")*results.parValue("delta_rho131");
  double r1m1_31_1_single_error_test = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho031")*results.parError("delta_rho031")+results.parValue("delta_rho131")*results.parValue("delta_rho131")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho131")*results.parError("delta_rho131")+2*results.parValue("delta_rho031")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho031")+2*results.parValue("delta_rho031")*results.parValue("delta_rho131")*results.covariance("rho1m11","rho1m10")+2*results.parValue("delta_rho031")*results.parValue("rho1m10")*results.covariance("rho1m11","delta_rho131")+2*results.parValue("rho1m11")*results.parValue("delta_rho131")*results.covariance("delta_rho031","rho1m10")+2*results.parValue("rho1m11")*results.parValue("rho1m10")*results.covariance("delta_rho031","delta_rho131")+2*results.parValue("rho1m10")*results.parValue("delta_rho131")*results.covariance("rho1m10","delta_rho131"));
  outfile << "r1m1_31_1_single_test" << "\t" << r1m1_31_1_single_value_test << "\t" << r1m1_31_1_single_error_test << "\n";
  
  double r1m1_3m1_1_single_value_test = results.parValue("rho1m11")*results.parValue("delta_rho03m1")+results.parValue("rho1m10")*results.parValue("delta_rho13m1");
  double r1m1_3m1_1_single_error_test = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho13m1")*results.parError("delta_rho13m1")+2*results.parValue("delta_rho03m1")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho03m1")+2*results.parValue("delta_rho03m1")*results.parValue("delta_rho13m1")*results.covariance("rho1m11","rho1m10")+2*results.parValue("delta_rho03m1")*results.parValue("rho1m10")*results.covariance("rho1m11","delta_rho13m1")+2*results.parValue("rho1m11")*results.parValue("delta_rho13m1")*results.covariance("delta_rho03m1","rho1m10")+2*results.parValue("rho1m11")*results.parValue("rho1m10")*results.covariance("delta_rho03m1","delta_rho13m1")+2*results.parValue("rho1m10")*results.parValue("delta_rho13m1")*results.covariance("rho1m10","delta_rho13m1"));
  outfile << "r1m1_3m1_1_single_test" << "\t" << r1m1_3m1_1_single_value_test << "\t" << r1m1_3m1_1_single_error_test << "\n";
  
  double r10_33_1_single_value_test = results.parValue("rho101")*(0.5-results.parValue("delta_rho011"))+results.parValue("rho100")*results.parValue("delta_rho133");
  double r10_33_1_single_error_test = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho011")*results.parError("delta_rho011")+results.parValue("delta_rho133")*results.parValue("delta_rho133")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho133")*results.parError("delta_rho133")-2*results.parValue("rho101")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho101","delta_rho011")+2*(0.5-results.parValue("delta_rho011"))*results.parValue("delta_rho133")*results.covariance("rho101","rho100")+2*results.parValue("rho101")*results.parValue("rho100")*results.covariance("rho101","delta_rho133")-2*results.parValue("rho101")*results.parValue("delta_rho133")*results.covariance("delta_rho011","rho100")-2*results.parValue("rho101")*results.parValue("rho100")*results.covariance("delta_rho011","delta_rho133")+2*results.parValue("rho100")*results.parValue("delta_rho133")*results.covariance("rho100","delta_rho133"));
  outfile << "r10_33_1_single_test" << "\t" << r10_33_1_single_value_test << "\t" << r10_33_1_single_error_test << "\n";

  double r10_31_1_single_value_test = results.parValue("rho101")*results.parValue("delta_rho031")+results.parValue("rho100")*results.parValue("delta_rho131");
  double r10_31_1_single_error_test = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho031")*results.parError("delta_rho031")+results.parValue("delta_rho131")*results.parValue("delta_rho131")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho131")*results.parError("delta_rho131")+2*results.parValue("delta_rho031")*results.parValue("delta_rho131")*results.covariance("rho101","rho100")+2*results.parValue("delta_rho031")*results.parValue("rho100")*results.covariance("rho101","delta_rho131")+2*results.parValue("rho101")*results.parValue("delta_rho131")*results.covariance("delta_rho031","rho100")+2*results.parValue("rho101")*results.parValue("rho100")*results.covariance("delta_rho031","delta_rho131")+2*results.parValue("delta_rho031")*results.parValue("rho101")*results.covariance("rho101","delta_rho031")+2*results.parValue("rho100")*results.parValue("delta_rho131")*results.covariance("rho100","delta_rho131"));
  outfile << "r10_31_1_single_test" << "\t" << r10_31_1_single_value_test << "\t" << r10_31_1_single_error_test << "\n";
  
  double r10_11_1_single_value_test = results.parValue("rho101")*results.parValue("delta_rho011")+results.parValue("rho100")*results.parValue("delta_rho111");
  double r10_11_1_single_error_test = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho011")*results.parError("delta_rho011")+results.parValue("delta_rho111")*results.parValue("delta_rho111")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho111")*results.parError("delta_rho111")+2*results.parValue("delta_rho011")*results.parValue("delta_rho111")*results.covariance("rho101","rho100")+2*results.parValue("delta_rho011")*results.parValue("rho100")*results.covariance("rho101","delta_rho111")+2*results.parValue("rho101")*results.parValue("delta_rho111")*results.covariance("delta_rho011","rho100")+2*results.parValue("rho101")*results.parValue("rho100")*results.covariance("delta_rho011","delta_rho111")+2*results.parValue("delta_rho011")*results.parValue("rho101")*results.covariance("rho101","delta_rho011")+2*results.parValue("rho100")*results.parValue("delta_rho111")*results.covariance("rho100","delta_rho111"));
  outfile << "r10_11_1_single_test" << "\t" << r10_11_1_single_value_test << "\t" << r10_11_1_single_error_test << "\n";
  
  double r10_3m1_1_single_value_test = results.parValue("rho101")*results.parValue("delta_rho03m1")+results.parValue("rho100")*results.parValue("delta_rho13m1");
  double r10_3m1_1_single_error_test = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+results.parValue("delta_rho13m1")*results.parValue("delta_rho13m1")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho13m1")*results.parError("delta_rho13m1")+2*results.parValue("delta_rho03m1")*results.parValue("delta_rho13m1")*results.covariance("rho101","rho100")+2*results.parValue("delta_rho03m1")*results.parValue("rho100")*results.covariance("rho101","delta_rho13m1")+2*results.parValue("rho101")*results.parValue("delta_rho13m1")*results.covariance("delta_rho03m1","rho100")+2*results.parValue("rho101")*results.parValue("rho100")*results.covariance("delta_rho03m1","delta_rho13m1")+2*results.parValue("delta_rho03m1")*results.parValue("rho101")*results.covariance("rho101","delta_rho03m1")+2*results.parValue("rho100")*results.parValue("delta_rho13m1")*results.covariance("rho100","delta_rho13m1"));
  outfile << "r10_3m1_1_single_test" << "\t" << r10_3m1_1_single_value_test << "\t" << r10_3m1_1_single_error_test << "\n";
 
  //factorized fitmodel 2
 //factorized fitmodel 2
double r00_33_1_single_value_test2 = results.parValue("rho001")*(0.5-results.parValue("delta_rho011")) ;
double r00_33_1_single_error_test2 = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho001")*results.parError("rho001")+results.parValue("rho001")*results.parValue("rho001")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho001")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho001","delta_rho011"));
outfile << "r00_33_1_single_test2" << "\t" << r00_33_1_single_value_test2 << "\t" << r00_33_1_single_error_test2 << "\n";

double r00_31_1_single_value_test2 = results.parValue("rho001")*results.parValue("delta_rho031") ;
double r00_31_1_single_error_test2 = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho001")*results.parError("rho001")+results.parValue("rho001")*results.parValue("rho001")*results.parError("delta_rho031")*results.parError("delta_rho031")+2*results.parValue("rho001")*results.parValue("delta_rho031")*results.covariance("rho001","delta_rho031"));
outfile << "r00_31_1_single_test2" << "\t" << r00_31_1_single_value_test2 << "\t" << r00_31_1_single_error_test2 << "\n";

double r00_11_1_single_value_test2 = results.parValue("rho001")*results.parValue("delta_rho011");
double r00_11_1_single_error_test2 = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho001")*results.parError("rho001")+results.parValue("rho001")*results.parValue("rho001")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("rho001")*results.parValue("delta_rho011")*results.covariance("rho001","delta_rho011"));
outfile << "r00_11_1_single_test2" << "\t" << r00_11_1_single_value_test2 << "\t" << r00_11_1_single_error_test2 << "\n";

double r00_3m1_1_single_value_test2 = results.parValue("rho001")*results.parValue("delta_rho03m1") ;
double r00_3m1_1_single_error_test2 = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho001")*results.parError("rho001")+results.parValue("rho001")*results.parValue("rho001")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+2*results.parValue("rho001")*results.parValue("delta_rho03m1")*results.covariance("rho001","delta_rho03m1"));
outfile << "r00_3m1_1_single_test2" << "\t" << r00_3m1_1_single_value_test2 << "\t" << r00_3m1_1_single_error_test2 << "\n";

double r11_33_1_single_value_test2 = results.parValue("rho111")*(0.5-results.parValue("delta_rho011"));
double r11_33_1_single_error_test2 = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho111")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho111","delta_rho011"));
outfile << "r11_33_1_single_test2" << "\t" << r11_33_1_single_value_test2 << "\t" << r11_33_1_single_error_test2 << "\n";

double r11_31_1_single_value_test2 = results.parValue("rho111")*results.parValue("delta_rho031");
double r11_31_1_single_error_test2 = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho031")*results.parError("delta_rho031")+2*results.parValue("rho111")*results.parValue("delta_rho031")*results.covariance("rho111","delta_rho031"));
outfile << "r11_31_1_single_test2" << "\t" << r11_31_1_single_value_test2 << "\t" << r11_31_1_single_error_test2 << "\n";

double r11_11_1_single_value_test2 = results.parValue("rho111")*results.parValue("delta_rho011");
double r11_11_1_single_error_test2 = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("rho111")*results.parValue("delta_rho011")*results.covariance("rho111","delta_rho011"));
outfile << "r11_11_1_single_test2" << "\t" << r11_11_1_single_value_test2 << "\t" << r11_11_1_single_error_test2 << "\n";

double r11_3m1_1_single_value_test2 = results.parValue("rho111")*results.parValue("delta_rho03m1");
double r11_3m1_1_single_error_test2 = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho111")*results.parError("rho111")+results.parValue("rho111")*results.parValue("rho111")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+2*results.parValue("rho111")*results.parValue("delta_rho03m1")*results.covariance("rho111","delta_rho03m1"));
outfile << "r11_3m1_1_single_test2" << "\t" << r11_3m1_1_single_value_test2 << "\t" << r11_3m1_1_single_error_test2 << "\n";

double r1m1_33_1_single_value_test2 = results.parValue("rho1m11")*(0.5-results.parValue("delta_rho011"));
double r1m1_33_1_single_error_test2 = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho1m11")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho1m11","delta_rho011"));
outfile << "r1m1_33_1_single_test2" << "\t" << r1m1_33_1_single_value_test2 << "\t" << r1m1_33_1_single_error_test2 << "\n";

double r1m1_11_1_single_value_test2 = results.parValue("rho1m11")*results.parValue("delta_rho011");
double r1m1_11_1_single_error_test2 = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("delta_rho011")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho011"));
outfile << "r1m1_11_1_single_test2" << "\t" << r1m1_11_1_single_value_test2 << "\t" << r1m1_11_1_single_error_test2 << "\n";

double r1m1_31_1_single_value_test2 = results.parValue("rho1m11")*results.parValue("delta_rho031");
double r1m1_31_1_single_error_test2 = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho031")*results.parError("delta_rho031")+2*results.parValue("delta_rho031")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho031"));
outfile << "r1m1_31_1_single_test2" << "\t" << r1m1_31_1_single_value_test2 << "\t" << r1m1_31_1_single_error_test2 << "\n";

double r1m1_3m1_1_single_value_test2 = results.parValue("rho1m11")*results.parValue("delta_rho03m1");
double r1m1_3m1_1_single_error_test2 = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho1m11")*results.parError("rho1m11")+results.parValue("rho1m11")*results.parValue("rho1m11")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+2*results.parValue("delta_rho03m1")*results.parValue("rho1m11")*results.covariance("rho1m11","delta_rho03m1"));
outfile << "r1m1_3m1_1_single_test2" << "\t" << r1m1_3m1_1_single_value_test2 << "\t" << r1m1_3m1_1_single_error_test2 << "\n";

double r10_33_1_single_value_test2 = results.parValue("rho101")*(0.5-results.parValue("delta_rho011"));
double r10_33_1_single_error_test2 = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho101")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho101","delta_rho011"));
outfile << "r10_33_1_single_test2" << "\t" << r10_33_1_single_value_test2 << "\t" << r10_33_1_single_error_test2 << "\n";

double r10_31_1_single_value_test2 = results.parValue("rho101")*results.parValue("delta_rho031");
double r10_31_1_single_error_test2 = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho031")*results.parError("delta_rho031")+2*results.parValue("delta_rho031")*results.parValue("rho101")*results.covariance("rho101","delta_rho031"));
outfile << "r10_31_1_single_test2" << "\t" << r10_31_1_single_value_test2 << "\t" << r10_31_1_single_error_test2 << "\n";

double r10_11_1_single_value_test2 = results.parValue("rho101")*results.parValue("delta_rho011");
double r10_11_1_single_error_test2 = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("delta_rho011")*results.parValue("rho101")*results.covariance("rho101","delta_rho011"));
outfile << "r10_11_1_single_test2" << "\t" << r10_11_1_single_value_test2 << "\t" << r10_11_1_single_error_test2 << "\n";

double r10_3m1_1_single_value_test2 = results.parValue("rho101")*results.parValue("delta_rho03m1");
double r10_3m1_1_single_error_test2 = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho101")*results.parError("rho101")+results.parValue("rho101")*results.parValue("rho101")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+2*results.parValue("delta_rho03m1")*results.parValue("rho101")*results.covariance("rho101","delta_rho03m1"));
outfile << "r10_3m1_1_single_test2" << "\t" << r10_3m1_1_single_value_test2 << "\t" << r10_3m1_1_single_error_test2 << "\n";



//allpha =2

  double r1m1_33_2_single_value_test = results.parValue("rho1m12")*(0.5-results.parValue("delta_rho011"));
  double r1m1_33_2_single_error_test = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho1m12")*results.parError("rho1m12")+results.parValue("rho1m12")*results.parValue("rho1m12")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho1m12")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho1m12","delta_rho011"));
  outfile << "r1m1_33_2_single_test" << "\t" << r1m1_33_2_single_value_test << "\t" << r1m1_33_2_single_error_test << "\n";

  double r1m1_11_2_single_value_test = results.parValue("rho1m12")*results.parValue("delta_rho011");
  double r1m1_11_2_single_error_test = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho1m12")*results.parError("rho1m12")+results.parValue("rho1m12")*results.parValue("rho1m12")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("rho1m12")*results.parValue("delta_rho011")*results.covariance("rho1m12","delta_rho011"));
  outfile << "r1m1_11_2_single_test" << "\t" << r1m1_11_2_single_value_test << "\t" << r1m1_11_2_single_error_test << "\n";

  double r1m1_31_2_single_value_test = results.parValue("rho1m12")*results.parValue("delta_rho031");///results.parValue("rho1m10")*results.parValue("delta_rho231");
  double r1m1_31_2_single_error_test = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho1m12")*results.parError("rho1m12")+results.parValue("rho1m12")*results.parValue("rho1m12")*results.parError("delta_rho031")*results.parError("delta_rho031")+results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho231")*results.parError("delta_rho231")+2*results.parValue("delta_rho031")*results.parValue("delta_rho231")*results.covariance("rho1m12","rho1m10")+2*results.parValue("delta_rho031")*results.parValue("rho1m10")*results.covariance("rho1m12","delta_rho231")+2*results.parValue("rho1m12")*results.parValue("delta_rho231")*results.covariance("delta_rho031","rho1m10")+2*results.parValue("rho1m12")*results.parValue("rho1m10")*results.covariance("delta_rho031","delta_rho231"));
  outfile << "r1m1_31_2_single_test" << "\t" << r1m1_31_2_single_value_test << "\t" << r1m1_31_2_single_error_test << "\n";

  double r1m1_3m1_2_single_value_test = results.parValue("rho1m12")*results.parValue("delta_rho03m1");//;
  double r1m1_3m1_2_single_error_test = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho1m12")*results.parError("rho1m12")+results.parValue("rho1m12")*results.parValue("rho1m12")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")+2*results.parValue("delta_rho03m1")*results.parValue("delta_rho23m1")*results.covariance("rho1m12","rho1m10")+2*results.parValue("delta_rho03m1")*results.parValue("rho1m10")*results.covariance("rho1m12","delta_rho23m1")+2*results.parValue("rho1m12")*results.parValue("delta_rho23m1")*results.covariance("delta_rho03m1","rho1m10")+2*results.parValue("rho1m12")*results.parValue("rho1m10")*results.covariance("delta_rho03m1","delta_rho23m1"));
  outfile << "r1m1_3m1_2_single_test" << "\t" << r1m1_3m1_2_single_value_test << "\t" << r1m1_3m1_2_single_error_test << "\n";

  double r10_33_2_single_value_test = results.parValue("rho102")*(0.5-results.parValue("delta_rho011"));
  double r10_33_2_single_error_test = sqrt((0.5-results.parValue("delta_rho011"))*(0.5-results.parValue("delta_rho011"))*results.parError("rho102")*results.parError("rho102")+results.parValue("rho102")*results.parValue("rho102")*results.parError("delta_rho011")*results.parError("delta_rho011")-2*results.parValue("rho102")*(0.5-results.parValue("delta_rho011"))*results.covariance("rho102","delta_rho011"));
  outfile << "r10_33_2_single_test" << "\t" << r10_33_2_single_value_test << "\t" << r10_33_2_single_error_test << "\n";

  double r10_11_2_single_value_test = results.parValue("rho102")*results.parValue("delta_rho011");
  double r10_11_2_single_error_test = sqrt(results.parValue("delta_rho011")*results.parValue("delta_rho011")*results.parError("rho102")*results.parError("rho102")+results.parValue("rho102")*results.parValue("rho102")*results.parError("delta_rho011")*results.parError("delta_rho011")+2*results.parValue("rho102")*results.parValue("delta_rho011")*results.covariance("rho102","delta_rho011"));
  outfile << "r10_11_2_single_test" << "\t" << r10_11_2_single_value_test << "\t" << r10_11_2_single_error_test << "\n";

  double r10_31_2_single_value_test = results.parValue("rho102")*results.parValue("delta_rho031"); //+results.parValue("rho100")*results.parValue("delta_rho231")
  double r10_31_2_single_error_test = sqrt(results.parValue("delta_rho031")*results.parValue("delta_rho031")*results.parError("rho102")*results.parError("rho102")+results.parValue("rho102")*results.parValue("rho102")*results.parError("delta_rho031")*results.parError("delta_rho031")+results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho231")*results.parError("delta_rho231")+2*results.parValue("delta_rho031")*results.parValue("delta_rho231")*results.covariance("rho102","rho100")+2*results.parValue("delta_rho031")*results.parValue("rho100")*results.covariance("rho102","delta_rho231")+2*results.parValue("rho102")*results.parValue("delta_rho231")*results.covariance("delta_rho031","rho100")+2*results.parValue("rho102")*results.parValue("rho100")*results.covariance("delta_rho031","delta_rho231"));
  outfile << "r10_31_2_single_test" << "\t" << r10_31_2_single_value_test << "\t" << r10_31_2_single_error_test << "\n";

  double r10_3m1_2_single_value_test = results.parValue("rho102")*results.parValue("delta_rho03m1"); //+results.parValue("rho100")*results.parValue("delta_rho23m1")
  double r10_3m1_2_single_error_test = sqrt(results.parValue("delta_rho03m1")*results.parValue("delta_rho03m1")*results.parError("rho102")*results.parError("rho102")+results.parValue("rho102")*results.parValue("rho102")*results.parError("delta_rho03m1")*results.parError("delta_rho03m1")+results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")+2*results.parValue("delta_rho03m1")*results.parValue("delta_rho23m1")*results.covariance("rho102","rho100")+2*results.parValue("delta_rho03m1")*results.parValue("rho100")*results.covariance("rho102","delta_rho23m1")+2*results.parValue("rho102")*results.parValue("delta_rho23m1")*results.covariance("delta_rho03m1","rho100")+2*results.parValue("rho102")*results.parValue("rho100")*results.covariance("delta_rho03m1","delta_rho23m1"));
  outfile << "r10_3m1_2_single_test" << "\t" << r10_3m1_2_single_value_test << "\t" << r10_3m1_2_single_error_test << "\n";

  double rt00_31_2_single_value_test = results.parValue("rho000")*results.parValue("delta_rho231");
  double rt00_31_2_single_error_test = sqrt(results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho000")*results.parError("rho000")+results.parValue("rho000")*results.parValue("rho000")*results.parError("delta_rho231")*results.parError("delta_rho231")+2*results.parValue("delta_rho231")*results.parValue("rho000")*results.covariance("rho000","delta_rho231"));
  outfile << "rt00_31_2_single_test" << "\t" << rt00_31_2_single_value_test << "\t" << rt00_31_2_single_error_test << "\n";

  double rt00_3m1_2_single_value_test = results.parValue("rho000")*results.parValue("delta_rho23m1");
  double rt00_3m1_2_single_error_test = sqrt(results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho000")*results.parError("rho000")+results.parValue("rho000")*results.parValue("rho000")*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")+2*results.parValue("delta_rho23m1")*results.parValue("rho000")*results.covariance("rho000","delta_rho23m1"));
  outfile << "rt00_3m1_2_single_test" << "\t" << rt00_3m1_2_single_value_test << "\t" << rt00_3m1_2_single_error_test << "\n";

  double rt11_31_2_single_value_test = (1-results.parValue("rho000"))*0.5*results.parValue("delta_rho231");
  double rt11_31_2_single_error_test = sqrt(0.25*results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho231")*results.parError("delta_rho231")-0.5*results.parValue("delta_rho231")*(1-results.parValue("rho000"))*results.covariance("rho000","delta_rho231"));
  outfile << "rt11_31_2_single_test" << "\t" << rt11_31_2_single_value_test << "\t" << rt11_31_2_single_error_test << "\n";

  double rt11_3m1_2_single_value_test = (1-results.parValue("rho000"))*0.5*results.parValue("delta_rho23m1");
  double rt11_3m1_2_single_error_test = sqrt(0.25*results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho000")*results.parError("rho000")+0.25*(1-results.parValue("rho000"))*(1-results.parValue("rho000"))*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")-0.5*results.parValue("delta_rho23m1")*(1-results.parValue("rho000"))*results.covariance("rho000","delta_rho23m1"));
  outfile << "rt11_3m1_2_single_test" << "\t" << rt11_3m1_2_single_value_test << "\t" << rt11_3m1_2_single_error_test << "\n";

  double rt1m1_31_2_single_value_test = results.parValue("rho1m10")*results.parValue("delta_rho231");
  double rt1m1_31_2_single_error_test = sqrt(results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho231")*results.parError("delta_rho231")+2*results.parValue("delta_rho231")*results.parValue("rho1m10")*results.covariance("rho1m10","delta_rho231"));
  outfile << "rt1m1_31_2_single_test" << "\t" << rt1m1_31_2_single_value_test << "\t" << rt1m1_31_2_single_error_test << "\n";

  double rt1m1_3m1_2_single_value_test = results.parValue("rho1m10")*results.parValue("delta_rho23m1");
  double rt1m1_3m1_2_single_error_test = sqrt(results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho1m10")*results.parError("rho1m10")+results.parValue("rho1m10")*results.parValue("rho1m10")*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")+2*results.parValue("delta_rho23m1")*results.parValue("rho1m10")*results.covariance("rho1m10","delta_rho23m1"));
  outfile << "rt1m1_3m1_2_single_test" << "\t" << rt1m1_3m1_2_single_value_test << "\t" << rt1m1_3m1_2_single_error_test << "\n";

  double rt10_31_2_single_value_test = results.parValue("rho100")*results.parValue("delta_rho231");
  double rt10_31_2_single_error_test = sqrt(results.parValue("delta_rho231")*results.parValue("delta_rho231")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho231")*results.parError("delta_rho231")+2*results.parValue("delta_rho231")*results.parValue("rho100")*results.covariance("rho100","delta_rho231"));
  outfile << "rt10_31_2_single_test" << "\t" << rt10_31_2_single_value_test << "\t" << rt10_31_2_single_error_test << "\n";

  double rt10_3m1_2_single_value_test = results.parValue("rho100")*results.parValue("delta_rho23m1");
  double rt10_3m1_2_single_error_test = sqrt(results.parValue("delta_rho23m1")*results.parValue("delta_rho23m1")*results.parError("rho100")*results.parError("rho100")+results.parValue("rho100")*results.parValue("rho100")*results.parError("delta_rho23m1")*results.parError("delta_rho23m1")+2*results.parValue("delta_rho23m1")*results.parValue("rho100")*results.covariance("rho100","delta_rho23m1"));
  outfile << "rt10_3m1_2_single_test" << "\t" << rt10_3m1_2_single_value_test << "\t" << rt10_3m1_2_single_error_test << "\n";
  

  // natural and unnatural components

double r00_33_0_plus_r00_33_1_value = results.parValue("r00_33_0") + results.parValue("r00_33_1");
double r00_33_0_plus_r00_33_1_error = sqrt(results.parError("r00_33_0")*results.parError("r00_33_0") + results.parError("r00_33_1")*results.parError("r00_33_1") + 2*results.covariance("r00_33_0","r00_33_1"));
outfile << "r00_33_0_plus_r00_33_1" << "\t" << r00_33_0_plus_r00_33_1_value << "\t" << r00_33_0_plus_r00_33_1_error << "\n";

double r00_33_0_minus_r00_33_1_value = results.parValue("r00_33_0") - results.parValue("r00_33_1");
double r00_33_0_minus_r00_33_1_error = sqrt(results.parError("r00_33_0")*results.parError("r00_33_0") + results.parError("r00_33_1")*results.parError("r00_33_1") - 2*results.covariance("r00_33_0","r00_33_1"));
outfile << "r00_33_0_minus_r00_33_1" << "\t" << r00_33_0_minus_r00_33_1_value << "\t" << r00_33_0_minus_r00_33_1_error << "\n";

double r00_11_0_plus_r00_11_1_value = 0.5 - results.parValue("r00_33_0") - 2.0*results.parValue("r11_11_0") - 2.0*results.parValue("r11_33_0") + results.parValue("r00_11_1");
double r00_11_0_plus_r00_11_1_error = sqrt(results.parError("r00_33_0")*results.parError("r00_33_0") + 4*results.parError("r11_11_0")*results.parError("r11_11_0") + 4*results.parError("r11_33_0")*results.parError("r11_33_0") + results.parError("r00_11_1")*results.parError("r00_11_1") + 4*results.covariance("r00_33_0","r11_11_0") + 4*results.covariance("r00_33_0","r11_33_0") - 2*results.covariance("r00_33_0","r00_11_1") + 8*results.covariance("r11_11_0","r11_33_0") - 4*results.covariance("r11_11_0","r00_11_1") - 4*results.covariance("r11_33_0","r00_11_1"));
outfile << "r00_11_0_plus_r00_11_1" << "\t" << r00_11_0_plus_r00_11_1_value << "\t" << r00_11_0_plus_r00_11_1_error << "\n";

double r00_11_0_minus_r00_11_1_value = 0.5 - results.parValue("r00_33_0") - 2.0*results.parValue("r11_11_0") - 2.0*results.parValue("r11_33_0") - results.parValue("r00_11_1");
double r00_11_0_minus_r00_11_1_error = sqrt(results.parError("r00_33_0")*results.parError("r00_33_0") + 4*results.parError("r11_11_0")*results.parError("r11_11_0") + 4*results.parError("r11_33_0")*results.parError("r11_33_0") + results.parError("r00_11_1")*results.parError("r00_11_1") + 4*results.covariance("r00_33_0","r11_11_0") + 4*results.covariance("r00_33_0","r11_33_0") + 2*results.covariance("r00_33_0","r00_11_1") + 8*results.covariance("r11_11_0","r11_33_0") + 4*results.covariance("r11_11_0","r00_11_1") + 4*results.covariance("r11_33_0","r00_11_1"));
outfile << "r00_11_0_minus_r00_11_1" << "\t" << r00_11_0_minus_r00_11_1_value << "\t" << r00_11_0_minus_r00_11_1_error << "\n";

double r00_31_0_plus_r00_31_1_value = results.parValue("r00_31_0") + results.parValue("r00_31_1");
double r00_31_0_plus_r00_31_1_error = sqrt(results.parError("r00_31_0")*results.parError("r00_31_0") + results.parError("r00_31_1")*results.parError("r00_31_1") + 2*results.covariance("r00_31_0","r00_31_1"));
outfile << "r00_31_0_plus_r00_31_1" << "\t" << r00_31_0_plus_r00_31_1_value << "\t" << r00_31_0_plus_r00_31_1_error << "\n";

double r00_31_0_minus_r00_31_1_value = results.parValue("r00_31_0") - results.parValue("r00_31_1");
double r00_31_0_minus_r00_31_1_error = sqrt(results.parError("r00_31_0")*results.parError("r00_31_0") + results.parError("r00_31_1")*results.parError("r00_31_1") - 2*results.covariance("r00_31_0","r00_31_1"));
outfile << "r00_31_0_minus_r00_31_1" << "\t" << r00_31_0_minus_r00_31_1_value << "\t" << r00_31_0_minus_r00_31_1_error << "\n";

double r00_3m1_0_plus_r00_3m1_1_value = results.parValue("r00_3m1_0") + results.parValue("r00_3m1_1");
double r00_3m1_0_plus_r00_3m1_1_error = sqrt(results.parError("r00_3m1_0")*results.parError("r00_3m1_0") + results.parError("r00_3m1_1")*results.parError("r00_3m1_1") + 2*results.covariance("r00_3m1_0","r00_3m1_1"));
outfile << "r00_3m1_0_plus_r00_3m1_1" << "\t" << r00_3m1_0_plus_r00_3m1_1_value << "\t" << r00_3m1_0_plus_r00_3m1_1_error << "\n";

double r00_3m1_0_minus_r00_3m1_1_value = results.parValue("r00_3m1_0") - results.parValue("r00_3m1_1");
double r00_3m1_0_minus_r00_3m1_1_error = sqrt(results.parError("r00_3m1_0")*results.parError("r00_3m1_0") + results.parError("r00_3m1_1")*results.parError("r00_3m1_1") - 2*results.covariance("r00_3m1_0","r00_3m1_1"));
outfile << "r00_3m1_0_minus_r00_3m1_1" << "\t" << r00_3m1_0_minus_r00_3m1_1_value << "\t" << r00_3m1_0_minus_r00_3m1_1_error << "\n";

double r11_33_0_plus_r1m1_33_1_value = results.parValue("r11_33_0") + results.parValue("r1m1_33_1");
double r11_33_0_plus_r1m1_33_1_error = sqrt(results.parError("r11_33_0")*results.parError("r11_33_0") + results.parError("r1m1_33_1")*results.parError("r1m1_33_1") + 2*results.covariance("r11_33_0","r1m1_33_1"));
outfile << "r11_33_0_plus_r1m1_33_1" << "\t" << r11_33_0_plus_r1m1_33_1_value << "\t" << r11_33_0_plus_r1m1_33_1_error << "\n";

double r11_33_0_minus_r1m1_33_1_value = results.parValue("r11_33_0") - results.parValue("r1m1_33_1");
double r11_33_0_minus_r1m1_33_1_error = sqrt(results.parError("r11_33_0")*results.parError("r11_33_0") + results.parError("r1m1_33_1")*results.parError("r1m1_33_1") - 2*results.covariance("r11_33_0","r1m1_33_1"));
outfile << "r11_33_0_minus_r1m1_33_1" << "\t" << r11_33_0_minus_r1m1_33_1_value << "\t" << r11_33_0_minus_r1m1_33_1_error << "\n";

double r11_11_0_plus_r1m1_11_1_value = results.parValue("r11_11_0") + results.parValue("r1m1_11_1");
double r11_11_0_plus_r1m1_11_1_error = sqrt(results.parError("r11_11_0")*results.parError("r11_11_0") + results.parError("r1m1_11_1")*results.parError("r1m1_11_1") + 2*results.covariance("r11_11_0","r1m1_11_1"));
outfile << "r11_11_0_plus_r1m1_11_1" << "\t" << r11_11_0_plus_r1m1_11_1_value << "\t" << r11_11_0_plus_r1m1_11_1_error << "\n";

double r11_11_0_minus_r1m1_11_1_value = results.parValue("r11_11_0") - results.parValue("r1m1_11_1");
double r11_11_0_minus_r1m1_11_1_error = sqrt(results.parError("r11_11_0")*results.parError("r11_11_0") + results.parError("r1m1_11_1")*results.parError("r1m1_11_1") - 2*results.covariance("r11_11_0","r1m1_11_1"));
outfile << "r11_11_0_minus_r1m1_11_1" << "\t" << r11_11_0_minus_r1m1_11_1_value << "\t" << r11_11_0_minus_r1m1_11_1_error << "\n";

double r11_31_0_plus_r1m1_31_1_value = results.parValue("r11_31_0") + results.parValue("r1m1_31_1");
double r11_31_0_plus_r1m1_31_1_error = sqrt(results.parError("r11_31_0")*results.parError("r11_31_0") + results.parError("r1m1_31_1")*results.parError("r1m1_31_1") + 2*results.covariance("r11_31_0","r1m1_31_1"));
outfile << "r11_31_0_plus_r1m1_31_1" << "\t" << r11_31_0_plus_r1m1_31_1_value << "\t" << r11_31_0_plus_r1m1_31_1_error << "\n";

double r11_31_0_minus_r1m1_31_1_value = results.parValue("r11_31_0") - results.parValue("r1m1_31_1");
double r11_31_0_minus_r1m1_31_1_error = sqrt(results.parError("r11_31_0")*results.parError("r11_31_0") + results.parError("r1m1_31_1")*results.parError("r1m1_31_1") - 2*results.covariance("r11_31_0","r1m1_31_1"));
outfile << "r11_31_0_minus_r1m1_31_1" << "\t" << r11_31_0_minus_r1m1_31_1_value << "\t" << r11_31_0_minus_r1m1_31_1_error << "\n";

double r11_3m1_0_plus_r1m1_3m1_1_value = results.parValue("r11_3m1_0") + results.parValue("r1m1_3m1_1");
double r11_3m1_0_plus_r1m1_3m1_1_error = sqrt(results.parError("r11_3m1_0")*results.parError("r11_3m1_0") + results.parError("r1m1_3m1_1")*results.parError("r1m1_3m1_1") + 2*results.covariance("r11_3m1_0","r1m1_3m1_1"));
outfile << "r11_3m1_0_plus_r1m1_3m1_1" << "\t" << r11_3m1_0_plus_r1m1_3m1_1_value << "\t" << r11_3m1_0_plus_r1m1_3m1_1_error << "\n";

double r11_3m1_0_minus_r1m1_3m1_1_value = results.parValue("r11_3m1_0") - results.parValue("r1m1_3m1_1");
double r11_3m1_0_minus_r1m1_3m1_1_error = sqrt(results.parError("r11_3m1_0")*results.parError("r11_3m1_0") + results.parError("r1m1_3m1_1")*results.parError("r1m1_3m1_1") - 2*results.covariance("r11_3m1_0","r1m1_3m1_1"));
outfile << "r11_3m1_0_minus_r1m1_3m1_1" << "\t" << r11_3m1_0_minus_r1m1_3m1_1_value << "\t" << r11_3m1_0_minus_r1m1_3m1_1_error << "\n";

double r1m1_33_0_plus_r11_33_1_value = results.parValue("r1m1_33_0") + results.parValue("r11_33_1");
double r1m1_33_0_plus_r11_33_1_error = sqrt(results.parError("r1m1_33_0")*results.parError("r1m1_33_0") + results.parError("r11_33_1")*results.parError("r11_33_1") + 2*results.covariance("r1m1_33_0","r11_33_1"));
outfile << "r1m1_33_0_plus_r11_33_1" << "\t" << r1m1_33_0_plus_r11_33_1_value << "\t" << r1m1_33_0_plus_r11_33_1_error << "\n";

double r1m1_33_0_minus_r11_33_1_value = results.parValue("r1m1_33_0") - results.parValue("r11_33_1");
double r1m1_33_0_minus_r11_33_1_error = sqrt(results.parError("r1m1_33_0")*results.parError("r1m1_33_0") + results.parError("r11_33_1")*results.parError("r11_33_1") - 2*results.covariance("r1m1_33_0","r11_33_1"));
outfile << "r1m1_33_0_minus_r11_33_1" << "\t" << r1m1_33_0_minus_r11_33_1_value << "\t" << r1m1_33_0_minus_r11_33_1_error << "\n";

double r1m1_11_0_plus_r11_11_1_value = results.parValue("r1m1_11_0") + results.parValue("r11_11_1");
double r1m1_11_0_plus_r11_11_1_error = sqrt(results.parError("r1m1_11_0")*results.parError("r1m1_11_0") + results.parError("r11_11_1")*results.parError("r11_11_1") + 2*results.covariance("r1m1_11_0","r11_11_1"));
outfile << "r1m1_11_0_plus_r11_11_1" << "\t" << r1m1_11_0_plus_r11_11_1_value << "\t" << r1m1_11_0_plus_r11_11_1_error << "\n";

double r1m1_11_0_minus_r11_11_1_value = results.parValue("r1m1_11_0") - results.parValue("r11_11_1");
double r1m1_11_0_minus_r11_11_1_error = sqrt(results.parError("r1m1_11_0")*results.parError("r1m1_11_0") + results.parError("r11_11_1")*results.parError("r11_11_1") - 2*results.covariance("r1m1_11_0","r11_11_1"));
outfile << "r1m1_11_0_minus_r11_11_1" << "\t" << r1m1_11_0_minus_r11_11_1_value << "\t" << r1m1_11_0_minus_r11_11_1_error << "\n";

double r1m1_31_0_plus_r11_31_1_value = results.parValue("r1m1_31_0") + results.parValue("r11_31_1");
double r1m1_31_0_plus_r11_31_1_error = sqrt(results.parError("r1m1_31_0")*results.parError("r1m1_31_0") + results.parError("r11_31_1")*results.parError("r11_31_1") + 2*results.covariance("r1m1_31_0","r11_31_1"));
outfile << "r1m1_31_0_plus_r11_31_1" << "\t" << r1m1_31_0_plus_r11_31_1_value << "\t" << r1m1_31_0_plus_r11_31_1_error << "\n";

double r1m1_31_0_minus_r11_31_1_value = results.parValue("r1m1_31_0") - results.parValue("r11_31_1");
double r1m1_31_0_minus_r11_31_1_error = sqrt(results.parError("r1m1_31_0")*results.parError("r1m1_31_0") + results.parError("r11_31_1")*results.parError("r11_31_1") - 2*results.covariance("r1m1_31_0","r11_31_1"));
outfile << "r1m1_31_0_minus_r11_31_1" << "\t" << r1m1_31_0_minus_r11_31_1_value << "\t" << r1m1_31_0_minus_r11_31_1_error << "\n";

double r1m1_3m1_0_plus_r11_3m1_1_value = results.parValue("r1m1_3m1_0") + results.parValue("r11_3m1_1");
double r1m1_3m1_0_plus_r11_3m1_1_error = sqrt(results.parError("r1m1_3m1_0")*results.parError("r1m1_3m1_0") + results.parError("r11_3m1_1")*results.parError("r11_3m1_1") + 2*results.covariance("r1m1_3m1_0","r11_3m1_1"));
outfile << "r1m1_3m1_0_plus_r11_3m1_1" << "\t" << r1m1_3m1_0_plus_r11_3m1_1_value << "\t" << r1m1_3m1_0_plus_r11_3m1_1_error << "\n";

double r1m1_3m1_0_minus_r11_3m1_1_value = results.parValue("r1m1_3m1_0") - results.parValue("r11_3m1_1");
double r1m1_3m1_0_minus_r11_3m1_1_error = sqrt(results.parError("r1m1_3m1_0")*results.parError("r1m1_3m1_0") + results.parError("r11_3m1_1")*results.parError("r11_3m1_1") - 2*results.covariance("r1m1_3m1_0","r11_3m1_1"));
outfile << "r1m1_3m1_0_minus_r11_3m1_1" << "\t" << r1m1_3m1_0_minus_r11_3m1_1_value << "\t" << r1m1_3m1_0_minus_r11_3m1_1_error << "\n";

double r10_33_0_plus_r10_33_1_value = results.parValue("r10_33_0") + results.parValue("r10_33_1");
double r10_33_0_plus_r10_33_1_error = sqrt(results.parError("r10_33_0")*results.parError("r10_33_0") + results.parError("r10_33_1")*results.parError("r10_33_1") + 2*results.covariance("r10_33_0","r10_33_1"));
outfile << "r10_33_0_plus_r10_33_1" << "\t" << r10_33_0_plus_r10_33_1_value << "\t" << r10_33_0_plus_r10_33_1_error << "\n";

double r10_33_0_minus_r10_33_1_value = results.parValue("r10_33_0") - results.parValue("r10_33_1");
double r10_33_0_minus_r10_33_1_error = sqrt(results.parError("r10_33_0")*results.parError("r10_33_0") + results.parError("r10_33_1")*results.parError("r10_33_1") - 2*results.covariance("r10_33_0","r10_33_1"));
outfile << "r10_33_0_minus_r10_33_1" << "\t" << r10_33_0_minus_r10_33_1_value << "\t" << r10_33_0_minus_r10_33_1_error << "\n";

double r10_11_0_plus_r10_11_1_value = results.parValue("r10_11_0") + results.parValue("r10_11_1");
double r10_11_0_plus_r10_11_1_error = sqrt(results.parError("r10_11_0")*results.parError("r10_11_0") + results.parError("r10_11_1")*results.parError("r10_11_1") + 2*results.covariance("r10_11_0","r10_11_1"));
outfile << "r10_11_0_plus_r10_11_1" << "\t" << r10_11_0_plus_r10_11_1_value << "\t" << r10_11_0_plus_r10_11_1_error << "\n";

double r10_11_0_minus_r10_11_1_value = results.parValue("r10_11_0") - results.parValue("r10_11_1");
double r10_11_0_minus_r10_11_1_error = sqrt(results.parError("r10_11_0")*results.parError("r10_11_0") + results.parError("r10_11_1")*results.parError("r10_11_1") - 2*results.covariance("r10_11_0","r10_11_1"));
outfile << "r10_11_0_minus_r10_11_1" << "\t" << r10_11_0_minus_r10_11_1_value << "\t" << r10_11_0_minus_r10_11_1_error << "\n";

double r10_31_0_plus_r10_31_1_value = results.parValue("r10_31_0") + results.parValue("r10_31_1");
double r10_31_0_plus_r10_31_1_error = sqrt(results.parError("r10_31_0")*results.parError("r10_31_0") + results.parError("r10_31_1")*results.parError("r10_31_1") + 2*results.covariance("r10_31_0","r10_31_1"));
outfile << "r10_31_0_plus_r10_31_1" << "\t" << r10_31_0_plus_r10_31_1_value << "\t" << r10_31_0_plus_r10_31_1_error << "\n";

double r10_31_0_minus_r10_31_1_value = results.parValue("r10_31_0") - results.parValue("r10_31_1");
double r10_31_0_minus_r10_31_1_error = sqrt(results.parError("r10_31_0")*results.parError("r10_31_0") + results.parError("r10_31_1")*results.parError("r10_31_1") - 2*results.covariance("r10_31_0","r10_31_1"));
outfile << "r10_31_0_minus_r10_31_1" << "\t" << r10_31_0_minus_r10_31_1_value << "\t" << r10_31_0_minus_r10_31_1_error << "\n";

double r10_3m1_0_plus_r10_3m1_1_value = results.parValue("r10_3m1_0") + results.parValue("r10_3m1_1");
double r10_3m1_0_plus_r10_3m1_1_error = sqrt(results.parError("r10_3m1_0")*results.parError("r10_3m1_0") + results.parError("r10_3m1_1")*results.parError("r10_3m1_1") + 2*results.covariance("r10_3m1_0","r10_3m1_1"));
outfile << "r10_3m1_0_plus_r10_3m1_1" << "\t" << r10_3m1_0_plus_r10_3m1_1_value << "\t" << r10_3m1_0_plus_r10_3m1_1_error << "\n";

double r10_3m1_0_minus_r10_3m1_1_value = results.parValue("r10_3m1_0") - results.parValue("r10_3m1_1");
double r10_3m1_0_minus_r10_3m1_1_error = sqrt(results.parError("r10_3m1_0")*results.parError("r10_3m1_0") + results.parError("r10_3m1_1")*results.parError("r10_3m1_1") - 2*results.covariance("r10_3m1_0","r10_3m1_1"));
outfile << "r10_3m1_0_minus_r10_3m1_1" << "\t" << r10_3m1_0_minus_r10_3m1_1_value << "\t" << r10_3m1_0_minus_r10_3m1_1_error << "\n";

double rt10_31_0_plus_rt10_31_1_value = results.parValue("rt10_31_0") + results.parValue("rt10_31_1");
double rt10_31_0_plus_rt10_31_1_error = sqrt(results.parError("rt10_31_0")*results.parError("rt10_31_0") + results.parError("rt10_31_1")*results.parError("rt10_31_1") + 2*results.covariance("rt10_31_0","rt10_31_1"));
outfile << "rt10_31_0_plus_rt10_31_1" << "\t" << rt10_31_0_plus_rt10_31_1_value << "\t" << rt10_31_0_plus_rt10_31_1_error << "\n";

double rt10_31_0_minus_rt10_31_1_value = results.parValue("rt10_31_0") - results.parValue("rt10_31_1");
double rt10_31_0_minus_rt10_31_1_error = sqrt(results.parError("rt10_31_0")*results.parError("rt10_31_0") + results.parError("rt10_31_1")*results.parError("rt10_31_1") - 2*results.covariance("rt10_31_0","rt10_31_1"));
outfile << "rt10_31_0_minus_rt10_31_1" << "\t" << rt10_31_0_minus_rt10_31_1_value << "\t" << rt10_31_0_minus_rt10_31_1_error << "\n";

double rt10_3m1_0_plus_rt10_3m1_1_value = results.parValue("rt10_3m1_0") + results.parValue("rt10_3m1_1");
double rt10_3m1_0_plus_rt10_3m1_1_error = sqrt(results.parError("rt10_3m1_0")*results.parError("rt10_3m1_0") + results.parError("rt10_3m1_1")*results.parError("rt10_3m1_1") + 2*results.covariance("rt10_3m1_0","rt10_3m1_1"));
outfile << "rt10_3m1_0_plus_rt10_3m1_1" << "\t" << rt10_3m1_0_plus_rt10_3m1_1_value << "\t" << rt10_3m1_0_plus_rt10_3m1_1_error << "\n";

double rt10_3m1_0_minus_rt10_3m1_1_value = results.parValue("rt10_3m1_0") - results.parValue("rt10_3m1_1");
double rt10_3m1_0_minus_rt10_3m1_1_error = sqrt(results.parError("rt10_3m1_0")*results.parError("rt10_3m1_0") + results.parError("rt10_3m1_1")*results.parError("rt10_3m1_1") - 2*results.covariance("rt10_3m1_0","rt10_3m1_1"));
outfile << "rt10_3m1_0_minus_rt10_3m1_1" << "\t" << rt10_3m1_0_minus_rt10_3m1_1_value << "\t" << rt10_3m1_0_minus_rt10_3m1_1_error << "\n";
  

  
  // Constraint 

  // double rhoN00 = 1/2*(results.parValue("rho000") - rho001_value);

  // double rhoU00 = 1/2*(results.parValue("rho000") + rho001_value);

  // double rhoN00 = 1/2*(results.parValue("rho110") + rhom111_value);

  // double rhoU00 = 1/2*(results.parValue("rho000") + rho001_value);
   
  double sigma_low = 2*(results.parValue("delta_rho111")+results.parValue("delta_rho133"));
  double sigma_low_err = sqrt(
    4*results.parError("delta_rho111")*results.parError("delta_rho111") +
    4*results.parError("delta_rho133")*results.parError("delta_rho133") +
    8*results.covariance("delta_rho111","delta_rho133")
);
  double sigma_up =(results.parValue("rho001")+2*results.parValue("rho111"));
  double sigma_up_err = sqrt(
    results.parError("rho001")*results.parError("rho001") +
    4*results.parError("rho111")*results.parError("rho111") +
    4*results.covariance("rho001","rho111")
);
outfile << "Sigma_low" << "\t" <<  sigma_low << "\t" <<  sigma_low_err << "\t"  << "\n";
outfile << "Sigma_up" << "\t" <<  sigma_up_err << "\t" <<  sigma_up_err << "\t"  << "\n";


  // covariance matrix
  vector< vector< double > > covMatrix;
  covMatrix = results.errorMatrix();
  
  // double SigmaN = results.parValue(pars[3]) + results.parValue(pars[6]);
  // double SigmaN_err = covMatrix[5][5] + covMatrix[8][8] + 2*covMatrix[5][8];

  // double SigmaD = 0.5*(1 - results.parValue(pars[0])) + results.parValue(pars[2]);
  // double SigmaD_err = 0.5*0.5*covMatrix[2][2] + covMatrix[4][4] - 2*0.5*covMatrix[2][4];

  // double Sigma = SigmaN/SigmaD;
  // double Sigma_err = fabs(Sigma) * sqrt(SigmaN_err/SigmaN/SigmaN + SigmaD_err/SigmaD/SigmaD);
  // outfile << Sigma << "\t" << Sigma_err << "\t";

  

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

