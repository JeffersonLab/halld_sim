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
  enum { kProd_Ang = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kVecMass = 5, kVecPsMass = 6, kt = 7, kRecoilMass = 8, kProtonPsMass = 9, kRecoilPsMass = 10, kNumHists};

  VecPs3PiPlotGenerator( const FitResults& results, Option opt);
  VecPs3PiPlotGenerator( const FitResults& results );
  VecPs3PiPlotGenerator( );
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );
 
};

#endif

