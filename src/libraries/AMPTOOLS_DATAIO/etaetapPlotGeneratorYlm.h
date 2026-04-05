#if !(defined ETAETAPPLOTGENERATORYLM)
#define ETAETAPPLOTGENERATORYLM

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

//class FitResults;
class Kinematics;

class etaetapPlotGeneratorYlm : public PlotGenerator
{

public:

  etaetapPlotGeneratorYlm( const FitResults& results );

  enum {
    khm12 = 0, khm13, khm23, khm1, khm2, khm3, kdltz, cosT, phiAng, cosT_m23, cosT_phi,
    kNumHists
  };

private:

  double polAngle; // polarization angle in DEGREES
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

};

#endif
