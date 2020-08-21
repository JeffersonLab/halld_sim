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
  enum { kOmegaPiMass = 0, kCosTheta = 1, kPhi = 2, kCosThetaH = 3, kPhiH = 4, kProd_Ang = 5, kt = 6, kNumHists};

  OmegaPiPlotGenerator( const FitResults& results );
  OmegaPiPlotGenerator( );
    
  void projectEvent( Kinematics* kin );
 
private:
  
  void createHistograms( );
 
};

#endif

