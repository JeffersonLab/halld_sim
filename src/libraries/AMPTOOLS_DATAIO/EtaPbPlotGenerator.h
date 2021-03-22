#if !(defined ETAPBPLOTGENERATOR)
#define ETAPBPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class EtaPbPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kTheta = 0, kPhi, kt, kTheta_phi, kt_phi, kNumHists};
  
  EtaPbPlotGenerator( const FitResults& results );
  EtaPbPlotGenerator( );
  
  void projectEvent( Kinematics* kin );
  
private:
  
  void createHistograms( );
  
};

#endif
