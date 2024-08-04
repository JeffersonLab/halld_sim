#if !defined(COMPLEXCOEFF)
#define COMPLEXCOEFF

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
void GPUComplexCoeff_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble real, GDouble imag);
#endif


class ComplexCoeff : public UserAmplitude< ComplexCoeff >
{  
public:
  
  ComplexCoeff() : UserAmplitude< ComplexCoeff >() { }
  ComplexCoeff( const vector< string >& args );
  
  ~ComplexCoeff(){}
  
  string name() const { return "ComplexCoeff"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

#ifdef GPU_ACCELERATION
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO) const{
	  GPUComplexCoeff_exec(dimGrid, dimBlock, GPU_AMP_ARGS, m_value.real(), m_value.imag());
  };
  
  bool isGPUEnabled() const { return true; }
#endif

  void updatePar( const AmpParameter& par );

private:
	
  AmpParameter m_param1;
  AmpParameter m_param2;
  bool m_represReIm;

  complex< GDouble > m_value;

};

#endif
