#if !(defined ISOPSPLOTGENERATOR)
#define ISOPSPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class IsoPsPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kProd_Ang = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kIsoMass = 5, kIsoPsMass = 6, kt = 7, kRecoilMass = 8, kProtonPsMass = 9, kRecoilPsMass = 10, kNumHists};

  IsoPsPlotGenerator( const FitResults& results, Option opt);
  IsoPsPlotGenerator( const FitResults& results );
  IsoPsPlotGenerator( );
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );
 
};

#endif

