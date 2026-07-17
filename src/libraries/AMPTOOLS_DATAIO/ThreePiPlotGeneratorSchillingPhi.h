#if !(defined THREEPIPLOTGENERATORSCHILLINGPHI)
#define THREEPIPLOTGENERATORSCHILLINGPHI

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class ThreePiPlotGeneratorSchillingPhi : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kPhi1020Mass = 0, kY20, kReY21Psi, kReY22Psi, kReY21Phi, kReY22Phi, kY40,
         kCosTheta, kPhi, kCosThetaPhi, kCosThetaPsi, kBigPhi, kPsi,
         kt, kCosThetaPiMinus, kCosThetaPiPlus, kCosThetaPi0, kPhiPiPlus,
         kPhiPiMinus, kPhiPi0, kNumHists};
  
  ThreePiPlotGeneratorSchillingPhi( const FitResults& results );
  ThreePiPlotGeneratorSchillingPhi( );
    
  void projectEvent( Kinematics* kin, const string& reactionName );

private:
        
  void createHistograms();
  
};

#endif
