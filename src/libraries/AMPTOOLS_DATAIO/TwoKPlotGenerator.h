#if !(defined TWOKPLOTGENERATOR)
#define TWOKPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class TwoKPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k2PiMass = 0, kPiPCosTheta, kPhiPiPlus, kPhiPiMinus, kPhi, kphi, kPsi, kt, kKpLabTheta, kKpLabPhi, kpKp, kKmLabTheta, kKmLabPhi, kpKm, kPLabTheta, kPLabPhi, kpP, kNumHists};


  TwoKPlotGenerator( const FitResults& results );
  TwoKPlotGenerator( );

  void projectEvent( Kinematics* kin );
  
private:
        
  void createHistograms();
  
};

#endif
