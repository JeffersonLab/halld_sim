#if !(defined TWOLEPTONPLOTGENERATOR)
#define TWOLEPTONPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoLeptonPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k2LeptonMass = 0, kLeptonPCosTheta, kPhiLeptonPlus, kPhiLeptonMinus, kPhi, kphi, kPsi, kt, kNumHists};
  
  TwoLeptonPlotGenerator( const FitResults& results );
  TwoLeptonPlotGenerator( );

  void projectEvent( Kinematics* kin );
  
private:
        
  void createHistograms();
  
};

#endif
