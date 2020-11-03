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
  enum { kOmegaPiMass = 0, kOmegaMass = 1, kCosTheta = 2, kPhi = 3, kCosThetaH = 4, kPhiH = 5, kProd_Ang = 6, kt = 7, kNumHists};

  OmegaPiPlotGenerator( const FitResults& results );
  OmegaPiPlotGenerator( );
    
  void projectEvent( Kinematics* kin );
 
private:
  
  void createHistograms( );
 
};

#endif
