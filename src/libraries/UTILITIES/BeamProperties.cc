/*
 *  BeamProperties.cc
 *
 *  Contains histograms for beam properties to be used in event generation and fitting.  Source
 *  of beam properties is from CombremsGeneration, external ROOT file or CCDB (to be implemented). 
 *
 *  Created by Justin Stevens on 12/29/17
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"

#include "BeamProperties.h"
#include "CobremsGeneration.hh"

using namespace ccdb;
using namespace std;

BeamProperties::BeamProperties( TString configFile ) {

	// check if histograms already exist before re-creating
	gDirectory->cd("/");
	fluxVsEgamma = (TH1D*)gDirectory->Get("BeamProperties_FluxVsEgamma");
	polFracVsEgamma = (TH1D*)gDirectory->Get("BeamProperties_PolFracVsEgamma");
	if(!fluxVsEgamma || !polFracVsEgamma) 
		createHistograms(configFile);
}

void BeamProperties::createHistograms( TString configFile ) {

	// First parse configuration file
	mConfigFile = configFile;
	bool isParsed = parseConfig();
	if(!isParsed) exit(1);

	// Fill flux histogram based on config file
	if(mIsCCDBFlux) 
		fillFluxFromCCDB();
	else if(mIsROOTFlux) 
		fillFluxFromROOT();
	else {  // default to CobremsGeneration for flux and polarization
		generateCobrems(); 
		if(mIsPolFixed) fillPolFixed();	// allow user to override with fixed polarization, if desired
		return;
	}
	
	// Fill polarization histogram based on config file
	if (mIsCCDBPol) 
		fillPolFromCCDB();
	else if(mIsROOTPol)
		fillPolFromROOT();
	else 
		fillPolFixed();	
	
	return;
}

bool BeamProperties::parseConfig(){

	cout<<endl<<"BeamProperties: Parsing config file "<<mConfigFile.Data()<<endl;

	// start assuming parameters for CobremsGeneration
	mIsCCDBFlux = false;
	mIsCCDBPol = false;
	mIsROOTFlux = false;
	mIsROOTPol = false;
	mIsPolFixed = false;

	ifstream inputFile;
	inputFile.open(mConfigFile.Data());
	if (!inputFile.is_open()){
		cout << "BeamProperties ERROR:  Could not open configuration file: " << mConfigFile << endl;
		return false;
	}

	while(inputFile) { 
		
		string line;
		getline(inputFile,line);
		
		// parse the line into words (from AmpTools ConfigFileParser)
		vector<string> words;
		string word("");
		for (unsigned int j = 0; j < line.size(); j++){
			if (!isspace(line[j])){
				word += line[j];
				if ((j == (line.size()-1))&&(!word.empty())){
					words.push_back(word);
					word = "";
				}
			}
			else if (!word.empty()){
				words.push_back(word);
				word = "";
			}
		}	

		// skip blank or incomplete lines and comments
		if(words.size() < 2 || words[0][0] == '#') 
			continue;

		// 1st is variable name and 2nd word is value
		if(words[0].find("CCDB") != std::string::npos) {
			mRunNumber = atoi(words[1].data());
		}
		else if(words[0].find("PolarizationMagnitude") != std::string::npos) {
			mIsPolFixed = true;
			mBeamParametersMap.insert( std::pair<std::string,double>( words[0].data(), atof(words[1].data()) ));
		}
		else if(words[0].find("ROOTFlux") != std::string::npos) {
			if(words[1].find("ccdb") != std::string::npos) 
				mIsCCDBFlux = true;
			else {
				mIsROOTFlux = true;
				mBeamHistNameMap.insert( std::pair<std::string,std::string>( words[0].data(), words[1].data() ) );
			}
		}
		else if(words[0].find("ROOTPol") != std::string::npos) {
			if(words[1].find("ccdb") != std::string::npos) 
				mIsCCDBPol = true;
			else {
				mIsROOTPol = true;
				mBeamHistNameMap.insert( std::pair<std::string,std::string>( words[0].data(), words[1].data() ) );
			}
		}
		else 
			mBeamParametersMap.insert( std::pair<std::string,double>( words[0].data(), atof(words[1].data()) ));
	}
	inputFile.close();

	return true;
}

// create histograms for flux and polarization fraction using CobremsGeneration
void BeamProperties::generateCobrems(){

	// Set default parameters (in case not provided by config file)
	std::map<std::string,double> defaultParameters;
	defaultParameters.insert( std::pair<std::string,double>( "ElectronBeamEnergy", 11.6 ) );
	defaultParameters.insert( std::pair<std::string,double>( "CoherentPeakEnergy", 8.8 ) );
	defaultParameters.insert( std::pair<std::string,double>( "PhotonBeamLowEnergy", 3.0 ) );
	defaultParameters.insert( std::pair<std::string,double>( "PhotonBeamHighEnergy", 11.6 ) );
	defaultParameters.insert( std::pair<std::string,double>( "Emittance", 2.5e-9 ) );
	defaultParameters.insert( std::pair<std::string,double>( "RadiatorThickness", 50.e-6 ) );
	defaultParameters.insert( std::pair<std::string,double>( "CollimatorDiameter", 0.005 ) );
	defaultParameters.insert( std::pair<std::string,double>( "CollimatorDistance", 76.0 ) );
	
	// check that required parameters are provided
	const int nParameters = 8;
	string parameterNames[nParameters] = {"ElectronBeamEnergy", "CoherentPeakEnergy", "PhotonBeamLowEnergy", "PhotonBeamHighEnergy",  "Emittance", "RadiatorThickness", "CollimatorDiameter", "CollimatorDistance"};
	for(int i=0; i<nParameters; i++) {
		if(mBeamParametersMap.find(parameterNames[i].data()) == mBeamParametersMap.end()) {
			cout << "BeamProperties WANRNING:  generateCombrems parameter " << parameterNames[i] << " missing: using default = " << defaultParameters.at(parameterNames[i].data()) << endl;
			mBeamParametersMap.insert( std::pair<std::string,double>( parameterNames[i].data(), defaultParameters.at(parameterNames[i].data())) );
		}
	}

	// Set parameters from config file
        double Emax  = mBeamParametersMap.at("ElectronBeamEnergy");
        double Epeak = mBeamParametersMap.at("CoherentPeakEnergy");
	double Elow  = mBeamParametersMap.at("PhotonBeamLowEnergy");
	double Ehigh = mBeamParametersMap.at("PhotonBeamHighEnergy");

	// Create histograms
	int nBinsEgamma = 1000;
	fluxVsEgamma = new TH1D("BeamProperties_FluxVsEgamma", "Flux vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	polFracVsEgamma = new TH1D("BeamProperties_PolFracVsEgamma", "Polarization Fraction vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	TH1D *polFluxVsEgamma   = new TH1D("BeamProperties_PolFluxVsEgamma", "Polarized Flux vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	
	// Setup cobrems
	CobremsGeneration cobrems(Emax, Epeak);
	cobrems.setBeamEmittance(mBeamParametersMap.at("Emittance"));
	cobrems.setTargetThickness(mBeamParametersMap.at("RadiatorThickness"));
	cobrems.setCollimatorDiameter(mBeamParametersMap.at("CollimatorDiameter"));
	cobrems.setCollimatorDistance(mBeamParametersMap.at("CollimatorDistance"));
	cobrems.setCollimatedFlag(true);
	
	// Fill flux
	cobrems.setPolarizedFlag(0); // 0=total flux
	for(int i=1; i<=nBinsEgamma; i++){
		double x = fluxVsEgamma->GetBinCenter(i)/Emax;
		double y = 0;
		if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
		else y = cobrems.Rate_dNtdx(x);
		fluxVsEgamma->SetBinContent(i, y);
	}

	// Fill polarized flux
	cobrems.setPolarizedFlag(1); // 1=polarized flux
	for(int i=1;i<=nBinsEgamma; i++){
		double x = fluxVsEgamma->GetBinCenter(i)/Emax;
		double y = 0;
		if(Epeak<Elow) y = 0.;
		else y = cobrems.Rate_dNcdx(x);
		polFluxVsEgamma->SetBinContent(i, y);
	}
	
	// Polarization fraction from ratio
	polFracVsEgamma->Divide(polFluxVsEgamma, fluxVsEgamma);

	return;
}

// load ROOT histograms for flux from external file
void BeamProperties::fillFluxFromROOT() {

	cout<<endl<<"BeamProperties: Using flux from input ROOT histogram file"<<endl;

	// open files and check that they exist
	TFile *fFlux;
	if(mBeamHistNameMap.count("ROOTFluxFile")) {
		fFlux = TFile::Open(mBeamHistNameMap.at("ROOTFluxFile").data());
	}
	else {
		cout << "BeamProperties ERROR:  ROOT flux file name not defined in configuration file" << endl;
		exit(1);
	}
	if(!fFlux->IsOpen()) {
		cout << "BeamProperties ERROR:  Could not open ROOT flux " << mBeamHistNameMap.at("ROOTFluxFile").data() << ", file doesn't exist" << endl;
		exit(1);
	}
	
	// open histograms and check that they exist
	if(mBeamHistNameMap.count("ROOTFluxName")) {
		fluxVsEgamma = (TH1D*)fFlux->Get(mBeamHistNameMap.at("ROOTFluxName").data())->Clone("BeamProperties_FluxVsEgamma");
	}
	else {
		cout << "BeamProperties ERROR:  ROOT flux histogram name not defined in configuration file" << endl;
		exit(1);
	}
	if(!fluxVsEgamma) {
		cout << "BeamProperties ERROR:  ROOT flux " << mBeamHistNameMap.at("ROOTFluxFile").data() << ", histogram doesn't exist" << endl;
		exit(1);
	}
	
	// set energy range for event generation
	if(!mBeamParametersMap.count("PhotonBeamLowEnergy") || !mBeamParametersMap.count("PhotonBeamHighEnergy")) {
		cout << "BeamProperties ERROR:  PhotonBeamLowEnergy or PhotonBeamHighEnergy not specified for event generation" << endl;
		exit(1);
	}
	fluxVsEgamma->GetXaxis()->SetRangeUser(mBeamParametersMap.at("PhotonBeamLowEnergy"), mBeamParametersMap.at("PhotonBeamHighEnergy")); 

	// keep in memory after file is closed
	fluxVsEgamma->SetDirectory(gROOT);
	fFlux->Close();
	
        return;
}

// load ROOT histograms for polarization fraction from external file
void BeamProperties::fillPolFromROOT() {

	cout<<endl<<"BeamProperties: Using polarization fraction from input ROOT histogram file"<<endl;

	// open files and check that they exist
	TFile *fPol;
	if(mBeamHistNameMap.count("ROOTPolFile")) {
		fPol = TFile::Open(mBeamHistNameMap.at("ROOTPolFile").data());
	}
	else {
		cout << "BeamProperties ERROR:  ROOT polarization file name not defined in configuration file" << endl;
		exit(1);
	}
	if(!fPol->IsOpen()) {
		cout << "BeamProperties ERROR:  Could not open ROOT polarization " << mBeamHistNameMap.at("ROOTPolFile").data() << ", file doesn't exist" << endl;
		exit(1);
	}
	
	// open histograms and check that they exist
	if(mBeamHistNameMap.count("ROOTPolName")) {
		polFracVsEgamma = (TH1D*)fPol->Get(mBeamHistNameMap.at("ROOTPolName").data())->Clone("BeamProperties_PolFracVsEgamma");
	}
	else {
		cout << "BeamProperties ERROR:  ROOT polarization histogram name not defined in configuration file" << endl;
		exit(1);
	}
	if(!polFracVsEgamma) {
		cout << "BeamProperties ERROR:  ROOT polarization " << mBeamHistNameMap.at("ROOTPolFile").data() << " histogram doesn't exist" << endl;
		exit(1);
	}

	// check contents of polarization histogram, and remove problematic bins with polarization > 100%
	for(int i=0; i<polFracVsEgamma->GetNbinsX(); i++) 
		if(fabs(polFracVsEgamma->GetBinContent(i)) > 1.0) 
			polFracVsEgamma->SetBinContent(i, 0);

	// keep in memory after file is closed
	polFracVsEgamma->SetDirectory(gROOT);
	fPol->Close();

        return;
}

double BeamProperties::PSAcceptance(double Egamma, double norm, double min, double max) {
	
	if( (Egamma > 2*min) && (Egamma < min + max)){
		return norm*(1-2*min/Egamma);
	} 
	else if( Egamma >= min + max){
		return norm*(2*max/Egamma - 1);
	}
	
	return 0.;
}


// create histograms for flux from CCDB for specified run number  and polarization fraction
void BeamProperties::fillFluxFromCCDB() {

	cout<<endl<<"BeamProperties: Using flux from CCDB run "<<mRunNumber<<endl;

	// Parse to get run number
	string ccdb_home(getenv("JANA_CALIB_URL"));
	string variation(getenv("JANA_CALIB_CONTEXT"));
	//cout<<ccdb_home.data()<<" "<<variation.data()<<endl;
	
	// Generate calibration class
	auto_ptr<ccdb::Calibration> calib(ccdb::CalibrationGenerator::CreateCalibration(ccdb_home, mRunNumber)); //, variation));
	
	// Get PS acceptance from CCDB
	vector< vector<double> > psAccept;
        calib->GetCalib(psAccept, "/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept");

	// find acceptance function parameters
	double PSnorm = psAccept[0][0];
	double PSmin = psAccept[0][1];
	double PSmax = psAccept[0][2];

	// Get tagger energy parameters from CCDB
	vector< double > photon_endpoint;
	vector< vector<double> > tagh_scaled_energy, tagm_scaled_energy;
	calib->GetCalib(photon_endpoint, "PHOTON_BEAM/endpoint_energy");
	calib->GetCalib(tagh_scaled_energy, "PHOTON_BEAM/hodoscope/scaled_energy_range");
	calib->GetCalib(tagm_scaled_energy, "PHOTON_BEAM/microscope/scaled_energy_range");

	// Setup custom histogram for filling flux from CCDB
	double Elow  = mBeamParametersMap.at("PhotonBeamLowEnergy"); // histogram min energy
	double Ehigh = mBeamParametersMap.at("PhotonBeamHighEnergy"); // histogram max energy
	int highest_tagh = 0;
	vector<double> Elows_tagh;
	for(int i=tagh_scaled_energy.size()-1; i>=0; i--) {
		double energy = tagh_scaled_energy[i][1] * photon_endpoint[0];
		if(energy > Elow && energy < Ehigh) { // keep only untagged flux for requested range
			Elows_tagh.push_back(energy);
			highest_tagh = i;
		}
	}
	Elows_tagh.push_back(tagh_scaled_energy[highest_tagh][2] * photon_endpoint[0]);	// add high energy edge for last counter
	sort(Elows_tagh.begin(), Elows_tagh.end());
      	fluxVsEgamma = new TH1D("BeamProperties_FluxVsEgamma", "Flux vs. E_{#gamma}", Elows_tagh.size()-1, &Elows_tagh[0]);

	// Get untagged flux from CCDB and fill histogram
	vector< vector<double> > taghflux, tagmflux;
        calib->GetCalib(taghflux, "PHOTON_BEAM/pair_spectrometer/lumi/tagh/untagged");
	calib->GetCalib(tagmflux, "PHOTON_BEAM/pair_spectrometer/lumi/tagm/untagged");
	for(uint i=0; i<taghflux.size(); i++) {
		double energy = 0.5*(tagh_scaled_energy[i][1]+tagh_scaled_energy[i][2]) * photon_endpoint[0];
		double accept = PSAcceptance(energy, PSnorm, PSmin, PSmax);
		double flux = taghflux[i][1]/accept;
		if(accept < 0.01) flux = 0; // remove very low acceptance regions to avoid fluctuations
		fluxVsEgamma->Fill(energy, flux);
	}
	
	return;
}

// PLACEHOLDER: create histograms for polarization fraction from CCDB for specified run number
void BeamProperties::fillPolFromCCDB() {
	
	cout<<endl<<"BeamProperties WARNING:  Polarization from CCDB is not currently implemented.  Event sample with be generated unpolarized (i.e. polarization=0)."<<endl<<endl;
	
	double polMagnitude = 0.0;
	polFracVsEgamma = new TH1D("BeamProperties_PolFracVsEgamma", "Polarization Fraction vs. E_{#gamma}", 1, 0., 13.);
	polFracVsEgamma->SetBinContent(1, polMagnitude);
	
        return;
}

// create single-bin histograms for polarization fraction at fixed value for all energies
void BeamProperties::fillPolFixed() {

	double polMagnitude = 0.0;
	if(mBeamParametersMap.count("PolarizationMagnitude")) 
		polMagnitude = mBeamParametersMap.at("PolarizationMagnitude");	
	
	cout<<endl<<"BeamProperties: Using fixed polarization = "<<polMagnitude<<endl;

	polFracVsEgamma = (TH1D*)gDirectory->Get("BeamProperties_PolFracVsEgamma");
	if(!polFracVsEgamma) polFracVsEgamma = new TH1D("BeamProperties_PolFracVsEgamma", "Polarization Fraction vs. E_{#gamma}", 1, 0., 13.);
	polFracVsEgamma->SetBinContent(1, polMagnitude);
	
        return;
}

double BeamProperties::GetPolAngle() {
	if(mBeamParametersMap.count("PolarizationAngle"))
		return mBeamParametersMap.at("PolarizationAngle");
	else 
		cout<<endl<<"BeamProperties WARNING: Polarization angle requeseted by generator/fitter but not found in beam configuration file.  Using PolAngle = 0 degrees by default"<<endl<<endl;

	return 0.;
}
