#if !(defined KHYPERONPLOTGENERATOR)
#define KHYPERONPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class KHyperonPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kHyperonMass = 0, kCosThetax, kCosThetay, kCosThetaz, kCosThetaHyp, kphiHyp, kPhiHyp, kNumHists};
  
  KHyperonPlotGenerator( const FitResults& results );
  KHyperonPlotGenerator( );

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
  
private:
        
  void createHistograms();

  map< string, double > m_reactionAngleMap;
};

#endif
