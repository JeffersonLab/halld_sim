#if !(defined VCPS3PIPLOTGENERATOR)
#define VECPS3PIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class VecPs3PiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kVecPs3PiMass = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kProd_Ang = 5, kt = 6, kRecoilMass = 7, kProtonPsMass = 8, kRecoilPsMass = 9, kLambda = 10, kDalitz = 11,  kNumHists};

  VecPs3PiPlotGenerator( const FitResults& results, Option opt);
  VecPs3PiPlotGenerator( const FitResults& results );
  VecPs3PiPlotGenerator( );
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );
 
};

#endif

