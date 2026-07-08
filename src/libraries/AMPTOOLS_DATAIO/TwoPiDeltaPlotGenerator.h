#if !(defined TWOPIDELTAPLOTGENERATOR)
#define TWOPIDELTAPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoPiDeltaPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum {
    k2PiMass = 0,
    kPPipMass,
    kPPimMass,
    kPi0PimMass,
    kPPipPimMass,
    kPPipPi0Mass,
    kPipPi0PimMass,
    kPPimPi0Mass,
    kPipPimMass,
    kPPi0Mass,
    kPipPi0Mass,
    kThetaDelta,
    kThetaPiPlus,
    kThetaPiMinus,
    kMomPiPlus,
    kMomPiMinus,
    klongMomPiPlus,
    klongMomPiMinus,
    klongMomPi0,
    klongMomProton,
    kMomPi0,
    kMomProton,
    kPhiPiPlus,
    kPhiPiMinus,
    kPhiProton,
    kphi_PiMinus_hel_rho,
    kphi_PiMinus_GJ_rho,
    kphi_PiPlus_hel_Delta,
    kphi_PiPlus_GJ_Delta,
    kphi_Proton_hel_Delta,
    kphi_Proton_GJ_Delta,
    kCosTheta_PiMinus_hel_rho,
    kCosTheta_PiMinus_GJ_rho,
    kCosTheta_PiPlus_hel_Delta,
    kCosTheta_PiPlus_GJ_Delta,
    kCosTheta_Proton_hel_Delta,
    kCosTheta_Proton_GJ_Delta,
    kPhi,
    kPsi,
    kt,
    kNumHists
  };
  
  TwoPiDeltaPlotGenerator( const FitResults& results );
  TwoPiDeltaPlotGenerator( );

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
  
private:
        
  void createHistograms();

  map< string, double > m_reactionAngleMap;
};

#endif
