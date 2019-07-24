//July 19th 2018, Based on DOI: 10.1016/0550-3213(84)90382-1
#if !defined(OMEGAPIANGAMP)
#define OMEGAPIANGAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"

using std::complex;
using namespace std;

class Kinematics;

class omegapiAngAmp : public UserAmplitude< omegapiAngAmp >
{

public:
  
  omegapiAngAmp() : UserAmplitude< omegapiAngAmp >() { }
  omegapiAngAmp( const vector< string >& args );
  
  string name() const { return "omegapiAngAmp"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

 void updatePar( const AmpParameter& par );

private:
  
  AmpParameter m_1p;
  AmpParameter w_1p;
  AmpParameter n_1p;
  AmpParameter phi0_1p;
  AmpParameter phip_1p;
  AmpParameter phim_1p;
  AmpParameter theta_1p;
  AmpParameter psi_1p;

  AmpParameter m_1m;
  AmpParameter w_1m;
  AmpParameter n_1m;
  AmpParameter phi0_1m;
  AmpParameter phip_1m;
  AmpParameter phim_1m;
  AmpParameter theta_1m;
  AmpParameter psi_1m;

  AmpParameter m_0m;
  AmpParameter w_0m;
  AmpParameter n_0m;
  AmpParameter phi0_0m;
  AmpParameter theta_0m;

  AmpParameter ds_ratio;

  double polAngle, polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;

};

#endif
