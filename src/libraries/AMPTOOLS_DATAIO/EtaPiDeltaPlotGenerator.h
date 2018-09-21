#if !(defined ETAPIDELTAPLOTGENERATOR)
#define ETAPIDELTAPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class EtaPiDeltaPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  //enum { kEtaPiMass = 0, kDeltaPPMass=0, kEtaCosTheta, kPhi, kt, kNumHists};
    enum{kEtaPiMass = 0, kDeltaPPMass, kEtaCosTheta, kPhi, kt, kEtaCosTheta_vs_Mass, kNumHists};
  EtaPiDeltaPlotGenerator( const FitResults& results );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
