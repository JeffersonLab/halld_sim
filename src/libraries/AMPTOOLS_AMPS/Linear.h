#if !defined(LINEAR)
#define LINEAR

#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPULinear_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     GDouble real_p0, GDouble real_p1,
                     GDouble imag_p0, GDouble imag_p1 );
#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class Linear : public UserAmplitude< Linear >
{
public:

  Linear() : UserAmplitude< Linear >() { }

  Linear( const vector< string >& args );

  ~Linear(){}

  string name() const { return "Linear"; }

  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

  enum UserVars { kMass = 0, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }

  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

  bool needsUserVarsOnly() const { return true; }
  bool areUserVarsStatic() const { return false; }

  void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

  bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:

  pair< string, string > m_daughters;
  AmpParameter m_real_p0;
  AmpParameter m_real_p1;
  AmpParameter m_imag_p0;

  GDouble m_imag_p1;

};

#endif
