#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

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

#include "AMPTOOLS_DATAIO/IsoPsPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/PiPiSWaveAMPK.h"
#include "AMPTOOLS_AMPS/Iso_ps_refl.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"



#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"



void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Uniform() );
  AmpToolsInterface::registerAmplitude( Iso_ps_refl() );
  AmpToolsInterface::registerAmplitude( PhaseOffset() );
  AmpToolsInterface::registerAmplitude( ComplexCoeff() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerAmplitude( OmegaDalitz() );
  AmpToolsInterface::registerAmplitude( PiPiSWaveAMPK() );

  
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
    cout << "\tisops_plotter <results file name> -o <output root file name> -t <output text file name>" << endl << endl;
    return 0;
  }

  bool showGui = false;
  string resultsName(argv[1]);
  string outName = "isops_histos.root";
  string outparsName = "isops__fitpars.txt";
  
  for (int i = 2; i < argc; i++){

    string arg(argv[i]);

    if (arg == "-o")  outName = argv[++i];
    if (arg == "-t")  outparsName = argv[++i];
    if (arg == "-g")  showGui = true;

    
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -o <file>\t output root file path" << endl;
      cout << "\t -t <file>\t output text file path" << endl;
      cout << "\t -g <file>\t show GUI" << endl;  
      exit(1);
    }
  }


    // ************************
    // parse the command line parameters
    // ************************

  cout << "Fit results file name    = " << resultsName << endl;
  cout << "Output root file name    = " << outName << endl;
  cout << "Output text file name    = " << outparsName << endl << endl;

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

	//	IsoPsPlotGenerator plotGen( results, PlotGenerator::kNoGenMC ); // optional can be omitted
	IsoPsPlotGenerator plotGen( results ); 
	cout << " Initialized ati and PlotGen" << endl;

 
    // *************************
    // Define Amplitudes and Coherent sums
    // *****************************

    vector<string> reflname = {"Uniform","PosRefl", "NegRefl"};
    vector<string> amphistname = {"Flat"};    
    vector<string> pipiIsobar_amps = {"pipiIso_0-S","pipiIso_1+P-","pipiIso_1+P0","pipiIso_1+P+"};
    vector<string> rhoIsobar_amps = {"rhoIso_0-P", "rhoIso_1+S-", "rhoIso_1+S0", "rhoIso_1+S+", "rhoIso_1+D-", "rhoIso_1+D0", "rhoIso_1+D+", "rhoIso_2+D--", "rhoIso_2+D-", "rhoIso_2+D0", "rhoIso_2+D+", "rhoIso_2+D++"}; 
    vector<string> f2Isobar_amps = {"f2Iso_1+P-","f2Iso_1+P0","f2Iso_1+P+","f2Iso_2+P--","f2Iso_2+P-","f2Iso_2+P0","f2Iso_2+P+","f2Iso_2+P++","f2Iso_2-S--","f2Iso_2-S-","f2Iso_2-S0","f2Iso_2-S+","f2Iso_2-S++","f2Iso_2-D--","f2Iso_2-D-","f2Iso_2-D0","f2Iso_2-D+","f2Iso_2-D++"};

    amphistname.insert(amphistname.end(), pipiIsobar_amps.begin(), pipiIsobar_amps.end());
    amphistname.insert(amphistname.end(), rhoIsobar_amps.begin(), rhoIsobar_amps.end());
    amphistname.insert(amphistname.end(), f2Isobar_amps.begin(), f2Isobar_amps.end());

    vector<string> ampsumname = {"pipiIso_0-S","pipiIso_1+P","rhoIso_0-P","rhoIso_1+S","rhoIso_1+D","rhoIso_2+D","f2Iso_1+P","f2Iso_2+P","f2Iso_2-S","f2Iso_2-D"};

    
    
    // ************************
    // set up an output ROOT file to store histograms
    // ************************

     TFile* plotfile = new TFile( outName.c_str(), "recreate");
     TH1::AddDirectory(kFALSE);
     TDirectory* plotfiledir = plotfile->mkdir("Contributions");  

     
     // *************************
     // Loop over different polarization types 
    // *************************
     size_t nReactions = results.reactionList().size();
     std::map<std::string, TH1*> summedHists;

     
  for (unsigned int polType = 0; polType < nReactions; polType++) {

    
    string reactionName = results.reactionList()[polType];

    //WRITING SEPARATE FILES FOR DIFF POLS    
    //  outName = reactionName + ".root";
    //TFile* plotfile = new TFile( outName.c_str(), "recreate");
    //TH1::AddDirectory(kFALSE);

    //    string dirname = "Contributions_"+reactionName;
    //    TDirectory* plotfiledir = plotfile->mkdir(dirname.c_str());  

  
    plotGen.enableReaction( reactionName );
    vector<string> sums = plotGen.uniqueSums();
    vector<string> amps = plotGen.uniqueAmplitudes();
    cout << "Reaction " << reactionName << " enabled with " << sums.size() << " sums and " << amps.size() << " amplitudes" << endl;

    string locampname;
  
    // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
    for (unsigned int irefl = 0; irefl <= reflname.size(); irefl++){

      // loop over desired amplitudes 
      for (unsigned int iamp = 0; iamp <= amphistname.size(); iamp++ ) {

	// turn on all ampltiudes by default 
	for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) plotGen.enableAmp( jamp );
      
	if (iamp < amphistname.size()){
	  // amplitude name
	  locampname = amphistname[iamp];	

	  // turn on all sums by default
	  for (unsigned int jsum = 0; jsum < sums.size(); jsum++) plotGen.enableSum(jsum);
	
	  //Arrange the negative and positive reflectivity sums 
	  //Negative reflectivity: "ImagNegSign" & "RealPosSign", as in the config file
	  //Positive reflectivity: "RealNegSign" & "ImagPosSign", as in the config file
	   	
	  if (irefl < 3) {
	    for (unsigned int i = 0; i < sums.size(); i++){
	      if( reflname[irefl] == "NegRefl" ){
		if(sums[i].find("RealNegSign") != std::string::npos || sums[i].find("ImagPosSign") != std::string::npos || sums[i].find("UniBG") != std::string::npos){
		  plotGen.disableSum(i);
		  //cout<<"disable sum "<<sums[i]<<"\n";
		}
	      }
	      if( reflname[irefl] == "PosRefl" ){
		if(sums[i].find("ImagNegSign") != std::string::npos || sums[i].find("RealPosSign") != std::string::npos || sums[i].find("UniBG") != std::string::npos) {
		  //cout<<"disable sum "<<sums[i]<<"\n";
		  plotGen.disableSum(i);
		}
	      }
	      if( reflname[irefl] == "Uniform" ){
		if(sums[i].find("ImagNegSign") != std::string::npos || sums[i].find("RealPosSign") != std::string::npos
		   || sums[i].find("RealNegSign") != std::string::npos || sums[i].find("ImagPosSign") != std::string::npos) {
		  //cout<<"disable sum "<<sums[i]<<"\n";
		  plotGen.disableSum(i);
		}
	      }
	    }	 
	  }
	  
	  //Pick up the right amplitude by excluding all nonmatching amplitude names
	  for (unsigned int jamp = 0; jamp < amps.size(); jamp++ ) 
	    if( amps[jamp].find(locampname.data()) == std::string::npos )  plotGen.disableAmp( jamp );

	} // close if (iamp < amphistname.size())


	
	// GENERATE HISTOGRAMS BY MEANS OF PLOT GENERATOR	
    
      cout << "Looping over input data" << endl;
      // loop over data, accMC, and genMC
      for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
		bool singleData =  irefl == reflname.size() && iamp == amphistname.size();
		bool singleFlatWave = (irefl == 0 && iamp > 0) || (irefl > 0 && iamp == 0);
		
		// if (iplot == PlotGenerator::kGenMC) continue; // no acceptance correction
		if ( iplot == PlotGenerator::kData && !singleData ) continue; // only plot data once
		if ( iplot == PlotGenerator::kBkgnd && !singleData ) continue; // only plot background once
		if ( iplot == PlotGenerator::kAccMC && singleFlatWave ) continue; // only plot Flat wave once
		if ( iplot == PlotGenerator::kGenMC && singleFlatWave ) continue; // only plot Flat wave once
		
		
	// loop over different variables
	for (unsigned int ivar  = 0; ivar  < IsoPsPlotGenerator::kNumHists; ivar++){
	  
	  // set reaction name as a prefix to the histogram name if all reactions are saved separately
	  // set a common name or no name as a prefix to the histogram name if all reactions are summed up 
	  //string histname = reactionName;
	  string histname = "";
	  
	  
	  if (ivar == IsoPsPlotGenerator::kProd_Ang)  histname += "Prod_Ang";
	  else if (ivar == IsoPsPlotGenerator::kCosTheta)  histname += "CosTheta_GJ";
	  else if (ivar == IsoPsPlotGenerator::kPhi)  histname += "Phi_GJ";
	  else if (ivar == IsoPsPlotGenerator::kCosThetaH)  histname += "CosTheta_HF";
	  else if (ivar == IsoPsPlotGenerator::kPhiH)  histname += "Phi_HF";
	  else if (ivar == IsoPsPlotGenerator::kIsoMass)  histname += "MIso";
	  else if (ivar == IsoPsPlotGenerator::kIsoPsMass)  histname += "MIsoPs";
	  else if (ivar == IsoPsPlotGenerator::kt)  histname += "minust";
	  else if (ivar == IsoPsPlotGenerator::kRecoilMass)  histname += "ProtonPiplusL_M";
	  else if (ivar == IsoPsPlotGenerator::kProtonPsMass)  histname += "ProtonPiminus_M";
	  else if (ivar == IsoPsPlotGenerator::kRecoilPsMass)  histname += "ProtonPiplusLPiminus_M";
	  else continue;	  

	  if (iplot == PlotGenerator::kData) histname += "_data";
	  if (iplot == PlotGenerator::kBkgnd) histname += "_bkgnd";
	  if (iplot == PlotGenerator::kAccMC) histname += "_fit";
	  if (iplot == PlotGenerator::kGenMC) histname += "_gen";
	  
	  if (irefl < reflname.size()){
	    // get name of sum for naming histogram
	    string sumName = reflname[irefl];
	    histname += "_";
	    histname += sumName;
	  }
	  //	  if (iamp > 0 && iamp < amphistname.size()) {
	  if (iamp < amphistname.size()) {
	    // get name of amp for naming histogram  
	    histname += "_";
	    histname += amphistname[iamp];
	  }
	  
	  Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
	  TH1* thist = hist->toRoot();
	  thist->SetName(histname.c_str());

	  // sum histograms across reactions
	  if (summedHists.find(histname) == summedHists.end())
	    summedHists[histname] = (TH1*) thist->Clone();
	  else 
	    summedHists[histname]->Add(thist);

	  //write separate histogram for each reaction
	  //if (iamp > 0 && iamp < amphistname.size())  plotfiledir->cd();
	  //else  plotfile->cd();
	  //thist->Write();
	  
	  delete thist;
	}
      }
    } // end of loop over amplitudes
   } // end of loop over sum configurations
  }// end of loop over 'reactions'

  

  //Now, write (nonzero) summed up histograms to the root file
  for (auto it = summedHists.begin(); it != summedHists.end(); it++) {

    std::string hsummed_name = it->first;
    TH1* hsummed = it->second;    

    if (hsummed->Integral() == 0.) continue;
    else hsummed->SetName(hsummed_name.c_str());
    
    if (hsummed_name.find("+") != std::string::npos || hsummed_name.find("-") != std::string::npos)  plotfiledir->cd();              
    else  plotfile->cd();    
    hsummed->Write();
  }

  
  plotfile->Close(); 



  // CALCULATE INTENSITY FRACTIONS AND PHASE DIFFERENCES
  // All 'reactions' are already turned on
  // The next calculations only need to happen one time
  
  // model parameters
  // cout << "Checking Parameters" << endl;
  
  // parameters to check
  vector< string > pars;  
  //  pars.push_back("dsratio");
  
  //Text file for writing the parameters
  ofstream outfile;
  outfile.open( outparsName );

  for(unsigned int i = 0; i<pars.size(); i++) {
    double parValue = results.parValue( pars[i] );
    double parError = results.parError( pars[i] );
    outfile << parValue << "\t" << parError << "\t" << endl;
  }


  
  //Define vector with full amplitudes
  vector<string> fullamps = plotGen.fullAmplitudes();

  //Define vectors for combining amplitudes with unique JPLMeps 
  const int nAmps = amphistname.size();
  vector<string> ampPosRefl[nAmps];
  vector<string> ampNegRefl[nAmps];
  vector< pair<string,string> > phaseDiffNames;

  //Define vectors for combining amplitudes with unique JPLeps 
  const int nSums = ampsumname.size();
  vector<string> ampsumPosRefl[nSums];
  vector<string> ampsumNegRefl[nSums];


  for(unsigned int i = 0; i < fullamps.size(); i++){    

    // Split by reflectivities and grab contributions for every JPLM amplitude 
    for(int iamp=0; iamp<nAmps; iamp++) {
	    string locampname = amphistname[iamp];
	    
	    if (fullamps[i].find("::" + locampname) == std::string::npos) continue;

	    if (fullamps[i].find("ImagNegSign") != std::string::npos || fullamps[i].find("RealPosSign") != std::string::npos)
	      ampNegRefl[iamp].push_back(fullamps[i]);

	    if (fullamps[i].find("ImagPosSign") != std::string::npos || fullamps[i].find("RealNegSign") != std::string::npos)
	      ampPosRefl[iamp].push_back(fullamps[i]);
    }

    // Split by reflectivities and grab contributions for every JPL sum 
    for(int isum=0; isum<nSums; isum++) {
	    string locsumname = ampsumname[isum];
	    
	    if (fullamps[i].find("::" + locsumname) == std::string::npos) continue;

	    if (fullamps[i].find("ImagNegSign") != std::string::npos || fullamps[i].find("RealPosSign") != std::string::npos)
	      ampsumNegRefl[isum].push_back(fullamps[i]);
	      
	    if (fullamps[i].find("ImagPosSign") != std::string::npos || fullamps[i].find("RealNegSign") != std::string::npos)
	      ampsumPosRefl[isum].push_back(fullamps[i]);
     }
    
    
    // second loop over amplitudes to get phase difference names
    for(unsigned int j = i+1; j < fullamps.size(); j++){

      // leave only the Spring2017_PARA_0::ImagPosSign and Spring2017_PARA_0::ImagNegSign coherent sums
      if (fullamps[i].find("Spring2017_PARA_0") == std::string::npos || fullamps[j].find("Spring2017_PARA_0") == std::string::npos) continue;      
      if (fullamps[i].find("UniBG") != std::string::npos || fullamps[j].find("UniBG") != std::string::npos) continue;      
      if (fullamps[i].find("Real") != std::string::npos || fullamps[j].find("Real") != std::string::npos) continue;      
      
      // only pair amplitudes in the same coherent sums
      if (fullamps[i].find("ImagNegSign") != std::string::npos && fullamps[j].find("ImagNegSign") == std::string::npos) continue;
      if (fullamps[i].find("ImagPosSign") != std::string::npos && fullamps[j].find("ImagPosSign") == std::string::npos) continue;
	    
      phaseDiffNames.push_back( std::make_pair(fullamps[i], fullamps[j]) );
      //      std::cout << "PAIRED:\t" << fullamps[i] << "\t" << fullamps[j] << endl; 
    }
  }


  outfile << "TOTAL EVENTS = " << results.intensity().first << " +- " << results.intensity().second << endl;
  
  cout<<"Computing intensities of amplitudes"<<endl;
  for (unsigned int i = 0; i < fullamps.size(); i++){
    vector<string> useamp;
    useamp.push_back(fullamps[i]);
    outfile << "FIT FRACTION " << fullamps[i] << " = " << results.intensity(useamp).first / results.intensity().first << " +- " << results.intensity(useamp).second / results.intensity().first <<  endl;
  }
  
  cout<<"Computing phase differences"<<endl;
  for(unsigned int i = 0; i < phaseDiffNames.size(); i++) {
	  pair <double, double> phaseDiff = results.phaseDiff( phaseDiffNames[i].first, phaseDiffNames[i].second );
	  outfile << "PHASE DIFF " << phaseDiffNames[i].first << " " << phaseDiffNames[i].second << " " << phaseDiff.first << " " << phaseDiff.second << endl;
  }

  cout<<"Computing intensities of coherent sums"<<endl;
  for(int i = 0; i < nAmps; i++){
    if(ampPosRefl[i].empty()) continue;
    outfile << "FIT FRACTION (coherent sum) PosRefl " << amphistname[i] << " = " << results.intensity(ampPosRefl[i]).first / results.intensity().first << " +- "
	    << results.intensity(ampPosRefl[i]).second / results.intensity().first << endl;
     outfile << "FIT FRACTION (coherent sum) NegRefl " << amphistname[i] << " = " << results.intensity(ampNegRefl[i]).first / results.intensity().first << " +- "
	     << results.intensity(ampNegRefl[i]).second / results.intensity().first << endl;
  }

  for(int i = 0; i < nSums; i++){
    if(ampsumPosRefl[i].empty()) continue;
    outfile << "FIT FRACTION (coherent sum) PosRefl " << ampsumname[i] << " = " << results.intensity(ampsumPosRefl[i]).first / results.intensity().first << " +- "
	    << results.intensity(ampsumPosRefl[i]).second / results.intensity().first << endl;
     outfile << "FIT FRACTION (coherent sum) NegRefl " << ampsumname[i] << " = " << results.intensity(ampsumNegRefl[i]).first / results.intensity().first << " +- "
	     << results.intensity(ampsumNegRefl[i]).second / results.intensity().first << endl;
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
