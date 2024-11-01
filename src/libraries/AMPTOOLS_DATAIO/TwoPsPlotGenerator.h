#if !(defined TWOPSPLOTGENERATOR)
#define TWOPSPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoPsPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  //enum { k2PiMass = 0, kPPipMass, kPPimMass, kPiPCosTheta, kThetaPiPlus, kThetaPiMinus, kThetaProton, kMomPiPlus, kMomPiMinus, kMomProton, kPhiPiPlus, kPhiPiMinus, kPhiProton, kPhi, kphi, kPsi, kt, kNumHists};
  enum {
    k2PsMass = 0,
    kLambdaKMass,
    kLambdaPiMass,
    kPiCosTheta,
    kPhiK,
    kPhiPi,
    kPhiLambda,
    kThetaK,
    kThetaPi,
    kThetaLambda,
    kMomK,
    kMomPi,
    kMomLambda,
    kPhi,
    kphi,
    kPsi,
    kt,
    kNumHists
  };

  TwoPsPlotGenerator( const FitResults& results );
  TwoPsPlotGenerator( );

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
  
private:
        
  void createHistograms();

  map< string, double > m_reactionAngleMap;
};

#endif
