#if !(defined ETAPLOTGENERATOR)
#define ETAPLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

//class FitResults;
class Kinematics;

class etapiPlotGenerator : public PlotGenerator
{

public:

  etapiPlotGenerator( const FitResults& results );

  enum {
    khm12 = 0, khm13, khm23, khm1, khm2, khm3, kdltz, cosT, phiAng, PhiT, cosT_m23, Omega, cosT_phi, cosT_Phi, cosT_lab, phiAng_lab, cosT_m23_lab, phi_m23_lab,
    kNumHists
  };

private:

  double polAngle; // polarization angle in DEGREES
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

};

#endif
