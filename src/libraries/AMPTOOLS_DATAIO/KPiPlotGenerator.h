#if !(defined KPIPLOTGENERATOR)
#define KPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class KPiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kKPiMass = 0, kLambKMass, kLambPiMass, kKCosTheta, kThetaK, kThetaPi, kThetaLamb, kMomK, kMomPi, kMomLamb, kPhiK, kPhiPi, kPhiLamb, kPhi, kphi, kPsi, kt, kNumHists};
  
  KPiPlotGenerator( const FitResults& results );
  KPiPlotGenerator( );

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
  
private:
        
  void createHistograms();

  map< string, double > m_reactionAngleMap;
};

#endif
