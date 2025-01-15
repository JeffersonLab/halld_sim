#if !defined(POLYNOMIAL)
#define POLYNOMIAL

#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>


using std::complex;
using namespace std;

class Kinematics;

#ifdef GPU_ACCELERATION
void GPUPolynomial_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO);

#endif


class Polynomial : public UserAmplitude< Polynomial >
{
public:

  Polynomial() : UserAmplitude< Polynomial >() { }

  Polynomial( const vector< string >& args );

  ~Polynomial(){}

  string name() const { return "Polynomial"; }

  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

#ifdef GPU_ACCELERATION
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const{
    GPUPolynomial_exec(dimGrid, dimBlock, GPU_AMP_ARGS);
  };

  bool isGPUEnabled() const { return true; }
#endif


private:

  pair< string, string > m_daughters;
  AmpParameter m_mag_p0;
  AmpParameter m_mag_p1;
  AmpParameter m_mag_p2;
  AmpParameter m_phase_p0;
  AmpParameter m_phase_p1;
  AmpParameter m_phase_p2;

};

#endif
