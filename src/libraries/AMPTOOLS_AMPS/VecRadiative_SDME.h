#if !defined(VECRADIATIVE_SDME)
#define VECRADIATIVE_SDME

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include "TFile.h"
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUVecRadiative_SDME_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                              GDouble rho000, GDouble rho100, GDouble rho1m10,
                              GDouble rho111, GDouble rho001, GDouble rho101,
                              GDouble rho1m11, GDouble rho102, GDouble rho1m12 );
#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class VecRadiative_SDME : public UserAmplitude< VecRadiative_SDME >
{
    
 public:
   
  VecRadiative_SDME() : UserAmplitude< VecRadiative_SDME >() { };
  VecRadiative_SDME( const vector< string >& args );
  
  enum UserVars { kPolFrac = 0, kCosTheta, kSinSqTheta, kSin2Theta,
                  kBigPhi, kPhi, kNumUserVars };
  
  unsigned int numUserVars() const { return kNumUserVars; }
  
  string name() const { return "VecRadiative_SDME"; }
    
  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

  // we can calcualte everything we need from userVars block so allow
  // the framework to purge the four-vectors
  bool needsUserVarsOnly() const { return true; }
  
#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
  bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:

  AmpParameter m_rho000;
  AmpParameter m_rho100;
  AmpParameter m_rho1m10;
  
  AmpParameter m_rho111;
  AmpParameter m_rho001;
  AmpParameter m_rho101;
  AmpParameter m_rho1m11;

  AmpParameter m_rho102;
  AmpParameter m_rho1m12;

  AmpParameter m_polAngle;

  GDouble m_polFraction;
  TH1D* m_polFrac_vs_E;
};

#endif

