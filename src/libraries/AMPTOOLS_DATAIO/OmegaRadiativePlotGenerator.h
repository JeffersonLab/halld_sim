#if !(defined OMEGARADIATIVEPLOTGENERATOR)
#define OMEGARADIATIVEPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class OmegaRadiativePlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kOmegaMass = 0, kY20, kReY21Psi, kReY22Psi, kReY21Phi, kReY22Phi, kY40,
         kCosTheta, kPhi, kCosThetaPhi, kCosThetaPsi, kBigPhi, kPsi,
         kt, kThetaLabPi0, kThetaLabGamma, kPThetaLabPi0, kPThetaLabGamma,
         kNumHists};

  OmegaRadiativePlotGenerator( const FitResults& results );
  OmegaRadiativePlotGenerator( );
    
  void projectEvent( Kinematics* kin, const string& reaction );
 
private:
  
  void createHistograms( );
 
};

#endif
