#if !(defined VCPSPLOTGENERATOR)
#define VECPSPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class VecPsPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kVecPsMass = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kProd_Ang = 5, kt = 6, kRecoilMass = 7, kProtonPsMass = 8, kRecoilPsMass = 9, kLambda = 10, kDalitz = 11,  kNumHists};

  VecPsPlotGenerator( const FitResults& results, Option opt);
  VecPsPlotGenerator( const FitResults& results );
  VecPsPlotGenerator( );
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );
 
};

#endif

