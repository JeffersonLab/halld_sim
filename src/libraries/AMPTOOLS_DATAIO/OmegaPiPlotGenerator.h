#if !(defined OMEGAPIPLOTGENERATOR)
#define OMEGAPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class OmegaPiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kOmegaPiMass = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kProd_Ang = 5, kt = 6, kRecoilMass = 7, kTwoPiMass = 8, kProtonPiMass = 9, kRecoilPiMass = 10, kLambda = 11, kDalitz = 12, kPhiDelta = 13, kCosThetaDelta = 14, kOmegaPiAngles = 15, kOmegaHAngles = 16, kDeltaAngles = 17, kBigLittlePhi = 18, kBigLittlePhiDelta = 19, kOmegaPipMass = 20, kOmega2PiMass = 21, kNumHists};

  OmegaPiPlotGenerator( const FitResults& results, Option opt);
  OmegaPiPlotGenerator( const FitResults& results );
  OmegaPiPlotGenerator( );
    
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

private:
  
  void createHistograms( );
 
};

#endif

