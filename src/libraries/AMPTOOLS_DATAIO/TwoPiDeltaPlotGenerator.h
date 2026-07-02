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
    kPhi,
    kPhi_GFJ,
    kphi_GFJ_rho, 
    kphi_HF_Delta,
    kphi_HF_rho,
    kphi_GFJ_Delta,
    kPsi,
    kt,
    kCosTheta_GFJ_Delta,
    kCosTheta_GFJ_rho,
    kCosTheta_HF_Delta,
    kCosTheta_HF_rho,
    kBeamassymetrie,
    kBeamassymetrie_Delta,

    
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
