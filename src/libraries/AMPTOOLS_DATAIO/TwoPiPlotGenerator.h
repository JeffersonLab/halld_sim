#if !(defined TWOPIPLOTGENERATOR)
#define TWOPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoPiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k2PiMass = 0, kPPipMass, kPPimMass, kPiPCosTheta, kThetaPiPlus, kThetaPiMinus, kThetaProton, kMomPiPlus, kMomPiMinus, kMomProton, kPhiPiPlus, kPhiPiMinus, kPhiProton, kPhi, kphi, kPsi, kt, kNumHists};
  
  TwoPiPlotGenerator( const FitResults& results );
  TwoPiPlotGenerator( );

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
  
private:
        
  void createHistograms();

  map< string, double > m_reactionAngleMap;
};

#endif
