#if !(defined ETAPI0PLOTGENERATOR)
#define ETAPI0PLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

class Kinematics;

class EtaPi0PlotGenerator : public PlotGenerator
{

public:

  EtaPi0PlotGenerator( const FitResults& results );

  enum {
    khm12 = 0, khm13, khm23, kdltz,cosT,phiAng, PhiT, cosT_m23, Omega,
    kNumHists
  };

private:

  double polAngle; // polarization angle in DEGREES
  void projectEvent( Kinematics* kin );

};

#endif
