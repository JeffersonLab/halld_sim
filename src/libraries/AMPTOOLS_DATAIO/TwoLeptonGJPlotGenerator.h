#if !(defined TWOLEPTONGJPLOTGENERATOR)
#define TWOLEPTONGJPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoLeptonGJPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k2LeptonMass = 0, kLeptonPCosTheta, kPhiLeptonPlus, kPhiLeptonMinus, kPhi, kphi, kPsi, kt, kNumHists};
  
  TwoLeptonGJPlotGenerator( const FitResults& results );
  TwoLeptonGJPlotGenerator( );

  void projectEvent( Kinematics* kin );
  
private:
        
  void createHistograms();
  
};

#endif
