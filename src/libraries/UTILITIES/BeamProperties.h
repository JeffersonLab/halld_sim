#if !defined(BEAMPROPERTIES)
#define BEAMPROPERTIES

/*
 *  BeamProperties.cc
 *
 *  Contains histograms for beam properties to be used in event generation and fitting.  Source
 *  of beam properties is from CombremsGeneration, external ROOT file or CCDB (to be implemented). 
 *
 *  Created by Justin Stevens on 12/29/17
 */

#include <string>
#include <map>

#include "TH1.h"

#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"
#include "TRandom3.h"
#include <TTimeStamp.h>

typedef vector< vector<double> > flux_t;

class BeamProperties {
  
public:
  
  BeamProperties( TString configFile);

  inline TH1D* GetFlux() { return fluxVsEgamma; };
  inline TH1D* GetPolFrac() { return polFracVsEgamma; };
  double GetPolAngle();

  // flux stored by channel
  flux_t taghflux;
  flux_t tagmflux;
	
private:

  void createHistograms( TString configFile );
  bool parseConfig();
  void generateCobrems();
  void fillFluxFromROOT();
  void fillPolFromROOT();
  void fillFluxFromCCDB();
  void fillTaggedFluxFromCCDB();
  void fillPolFromCCDB();
  void fillPolFixed();
  double PSAcceptance(double Egamma, double norm, double min, double max);

  TString mConfigFile;
  std::map<std::string,double> mBeamParametersMap;
  std::map<std::string,std::string> mBeamHistNameMap;

  bool mIsCCDBTaggedFlux;
  bool mIsCCDBFlux, mIsCCDBPol;
  bool mIsROOTFlux, mIsROOTPol;
  bool mIsPolFixed;
  int mRunNumber;

  TH1D *fluxVsEgamma;
  TH1D *polFracVsEgamma;

};

#endif
