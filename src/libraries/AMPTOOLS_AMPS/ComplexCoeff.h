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
void GPUComplexCoeff_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble param1, GDouble param2, bool represReIm);
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
	  GPUComplexCoeff_exec(dimGrid, dimBlock, GPU_AMP_ARGS, m_param1, m_param2, m_represReIm);
  };
  
  bool isGPUEnabled() const { return true; }
#endif

  void updatePar( const AmpParameter& par );

private:
	
  AmpParameter m_param1;
  AmpParameter m_param2;
  bool m_represReIm;

#ifndef GPU_ACCELERATION
  complex< GDouble > m_value;
#endif

};

#endif
